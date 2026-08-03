#include "solver_lumped.h"
#include <iomanip>

//--------------------------------------------------------------
solver_lumped::solver_lumped(string a_name, string a_folder)
{
	name = a_name;
	input_folder_path = a_folder;
}

solver_lumped::~solver_lumped(){}

//--------------------------------------------------------------
void solver_lumped::initialization(double time_period)
{
	// setting sizes
	number_of_nodes = nodes.size();
	number_of_edges = edges.size();
	number_of_master = boundary_model_node.size();

	// clearing master boundary indices
	boundary_indices.clear();

	// heart rate
	//heart_rate = hr; // from Charlton2019
	//time_period = 60./heart_rate;

	// setting the par variables, converting from SI to non-SI for favourable conditioning
	set_non_SI_parameters();

	double E = elastance(0.);
	for(int i=0; i<number_of_nodes; i++)
	{
		nodes[i]->pressure.clear();
		nodes[i]->pressure.push_back(nodes[i]->pressure_initial);
		nodes[i]->p = nodes[i]->pres_ini_non_SI;
		nodes[i]->y = nodes[i]->p/E;
		nodes[i]->RBC_fi0Dn = fi_init_RBC_lum;
		nodes[i]->HBsat_0Dn = init_HB_sat_lum;
		nodes[i]->PlasmaO2_0Dn = init_PlasmaO2_lum;
		//CO2 stuff
		nodes[i]->CO2_pla_n = fi_init_CO2_pla;
		nodes[i]->CO2_rbc_n = fi_init_CO2_rbc;
		nodes[i]->HCO3_pla_n = fi_init_HCO3_pla;
		nodes[i]->HCO3_rbc_n = fi_init_HCO3_rbc;
		nodes[i]->HbCO2_n = fi_init_HbCO2;
	}

	// building model
	build_system();

	number_of_elastance = 0;
	for(int i=0; i<number_of_edges; i++)
	{
		edges[i]->volume_flow_rate.clear();
		edges[i]->volume_flow_rate.push_back(edges[i]->volume_flow_rate_initial);
		edges[i]->vfr = edges[i]->vfr_ini_non_SI;
		if(edges[i]->type_code == 2)
		{
			number_of_elastance++;
			int si = edges[i]->node_index_start;
			int ei = edges[i]->node_index_end;
			E = elastance(0.,edges[i]->par_non_SI);
			nodes[si]->y = nodes[si]->p/E;
			nodes[ei]->y = nodes[ei]->p/E;
		}
	}

	// setting first time stamp
	time.clear();
	time.push_back(0.);

	//tissueO2 init
	tissueO2_save.clear();
	tissueO2_save.push_back(init_tissueO2);

	// setting Eigen vars
	int nm = number_of_edges + number_of_nodes + number_of_master + 2*number_of_elastance;
	A = MatrixXd::Zero(nm,nm);
	b = VectorXd::Zero(nm);


	// for myogenic control
	p_ave = new time_average();


	//for metabolic response
	Ct_ave = new time_average();

	// for CO2 response
	P_CO2_ave = new time_average();
}

//--------------------------------------------------------------
void solver_lumped::set_newton_size()
{
	number_of_moc = boundary_indices.size();

	// setting Eigen vars for nonlinear solvr
	int N = number_of_edges + number_of_nodes + 2*number_of_elastance + number_of_moc;
	Jac = MatrixXd::Zero(N,N);
	x = VectorXd::Zero(N);
	f = VectorXd::Zero(N);
}

//--------------------------------------------------------------
void solver_lumped::coefficients_newton(double t_act)
{
	// increasing time
	double dt = t_act - time.back();

	// sizes of nodes and edges
	int n=number_of_nodes, m=number_of_edges, l=number_of_elastance;

	// tracing the virtual nodes of elastance
	int i_elas=0;

	// edges
	for(int i=0; i<number_of_edges; i++)
	{
		int i1 = edges[i]->node_index_start;
		int i2 = edges[i]->node_index_end;

		double par = edges[i]->par_non_SI[0]*edges[i]->parameter_factor;

		if(edges[i]->type_code == 0) // resistor
		{
			Jac(i,m+i2) = 1.;
			Jac(i,m+i1) = -1.;
			Jac(i,i) = par; // R*Rf

			f(i) = x(m+i2) - x(m+i1) + par*x(i);
		}
		else if(edges[i]->type_code == 1) // capacitor
		{
			Jac(i,m+i2) = 1.;
			Jac(i,m+i1) = -1.;
			Jac(i,i) = dt/par; // dt/C

			f(i) = x(m+i2) - x(m+i1) + dt/par*x(i) - nodes[i2]->p + nodes[i1]->p;
		}
		else if(edges[i]->type_code == 2) // elastance
		{
			// actual elastance
			double E_act = elastance(t_act,edges[i]->par_non_SI);

			// basic equation for the edge
			Jac(i,m+n+i_elas+1) = 1.;
			Jac(i,m+n+i_elas) = -1.;
			Jac(i,i) = dt;
			f(i) = x(m+n+i_elas+1) - x(m+n+i_elas) + dt*x(i) - nodes[i2]->y + nodes[i1]->y;

			// equations for the virtual nodes
			Jac(n+m+i_elas,m+i1) = -1.;
			Jac(n+m+i_elas+1,m+i2) = -1.;
			Jac(n+m+i_elas,n+m+i_elas) = E_act;
			Jac(n+m+i_elas+1,n+m+i_elas+1) = E_act;
			f(n+m+i_elas) = E_act*x(n+m+i_elas) - x(m+i1);
			f(n+m+i_elas+1) = E_act*x(n+m+i_elas+1) - x(m+i2);
			i_elas+=2;
		}
		else if(edges[i]->type_code == 3) // inductor
		{
			Jac(i,m+i2) = 1.;
			Jac(i,m+i1) = -1.;
			Jac(i,i) = par/dt; // L/dt
			f(i) = x(m+i2) - x(m+i1) + par/dt * (x(i)-edges[i]->vfr);
		}
		else if(edges[i]->type_code == 4) // voltage source
		{
			Jac(i,m+i2) = 1.;
			Jac(i,m+i1) = -1.;
			f(i) = x(m+i2) - x(m+i1) - par;
		}
		else if(edges[i]->type_code == 5) // diode
		{
			Jac(i,m+i2) = 1.;
			Jac(i,m+i1) = -1.;

			if(x(m+i1)>x(m+i2)) // diode is open
			{
				edges[i]->is_open = true;
				Jac(i,i) = par; // R
				f(i) = x(m+i2) - x(m+i1) + par*x(i);
			}
			else // diode is closed
			{
				edges[i]->is_open = false;
				Jac(i,i) = 1.e10*par; // R*1.e10
				f(i) = x(m+i2) - x(m+i1) + 1.e10*par*x(i);
			}
		}
		else if(edges[i]->type_code == 6) // valve
		{
			Jac(i,m+i2) = 1.;
			Jac(i,m+i1) = -1.;
			Jac(i,i) = 2.*par*x(i); // 2*R*Q

			f(i) = x(m+i2) - x(m+i1) + par*x(i)*x(i); // dp = R*Q^2
		}
		else if(edges[i]->type_code == 7) // resistor for coronaries
		{
			double E_act = elastance(t_act);

			Jac(i,m+i2) = 1.;
			Jac(i,m+i1) = -1.;
			double R = par*(1. + beta_coronary*E_act/elastance_max_nom); // from Reymond2009
			Jac(i,i) = R; // R

			f(i) = x(m+i2) - x(m+i1) + R*x(i);
		}
		else if(edges[i]->type_code == 8) // capacitor
		{
			double E_act = elastance(t_act);

			Jac(i,m+i2) = 1.;
			Jac(i,m+i1) = -1.;
			double C = par*(1. - alpha_coronary*E_act/elastance_max_nom); // from Reymond2009
			Jac(i,i) = dt/C; // dt/C

			f(i) = x(m+i2) - x(m+i1) + dt/C*x(i) - nodes[i2]->p + nodes[i1]->p;
		}
		else if(edges[i]->type_code == 9) // current source
		{
			Jac(i,i) = 1.;
			f(i) = x(i) - par;
		}
		else if(edges[i]->type_code == 10) // piecevise constant resistor, models compression
		{

			//if(edges[i]->t1 > t_act){par = edges[i]->par_non_SI[0];}
			if (edges[i]->t1 < t_act && t_act < edges[i]->t2) 
			{
    			par = edges[i]->par_non_SI[1];
			}
			//else {par = edges[i]->par_non_SI[0];}

			Jac(i,m+i2) = 1.;
			Jac(i,m+i1) = -1.;
			Jac(i,i) = par; // R*Rf

			f(i) = x(m+i2) - x(m+i1) + par*x(i);
		}
	}

	// nodes
	for(int i=0; i<number_of_nodes; i++)
	{
		if(nodes[i]->is_ground == false) // intersections
		{
			f(m+i) = 0.;
			for(int j=0; j<nodes[i]->edge_in.size(); j++)
			{
				Jac(m+i,nodes[i]->edge_in[j]) = 1;
				f(m+i) += x(nodes[i]->edge_in[j]);
			}
			for(int j=0; j<nodes[i]->edge_out.size(); j++)
			{
				Jac(m+i,nodes[i]->edge_out[j]) = -1;
				f(m+i) -= x(nodes[i]->edge_out[j]);
			}
		}
		else // ground nodes, pi = p0[mmHg]
		{
			Jac(m+i,m+i) = 1;
			f(m+i) = x(m+i)-1.e5/mmHg_to_Pa;
		}
	}
}

//--------------------------------------------------------------
void solver_lumped::initialization_newton(double t_act)
{

	int i_elas=0;
	for(int i=0; i<number_of_edges; i++)
	{
		x(i) = edges[i]->vfr;
		if(edges[i]->type_code == 2) // elastance
		{
			int i1 = edges[i]->node_index_start;
			int i2 = edges[i]->node_index_end;
			x(number_of_edges+number_of_nodes+i_elas) = nodes[i1]->y;
			x(number_of_edges+number_of_nodes+i_elas+1) = nodes[i2]->y;
			i_elas+=2;
		}
	}
	for(int i=0; i<number_of_nodes; i++)
	{
		x(number_of_edges+i) = nodes[i]->p;
	}
}

//--------------------------------------------------------------
void solver_lumped::substitute_newton(double t_act)
{
	// saving time step
	time.push_back(t_act);

	// putting back the outputs
	vector<double> par{elastance_max_nom,elastance_min_nom};
	double E_act = elastance(time.back(),par);
	for(int i=0; i<number_of_nodes; i++)
	{
		nodes[i]->p = x(number_of_edges+i);
		nodes[i]->y = x(number_of_edges+i)/E_act;
		if(nodes[i]->do_save_memory)
		{
			nodes[i]->pressure.push_back(x(number_of_edges+i)*mmHg_to_Pa);
		}
	}

	for(int i=0; i<number_of_edges; i++)
	{
		edges[i]->vfr = x(i);
		if(edges[i]->do_save_memory)
		{
			edges[i]->volume_flow_rate.push_back(x(i)*1.e-6);
		}
		if(edges[i]->type_code == 2) // elastance
		{
			// rewriting the elastance nodes with actual E_act
			int si = edges[i]->node_index_start;
			int ei = edges[i]->node_index_end;
			E_act = elastance(time.back(),edges[i]->par_non_SI);
			nodes[si]->y = nodes[si]->p/E_act;
			nodes[ei]->y = nodes[ei]->p/E_act;
		}
	}


	// saving time averages of myogenic control
	if(do_myogenic)
	{
		double tn = time.back();
		//double vn = edges[0]->vfr;

		double vn = nodes[5]->p;
		//double vn = nodes[1]->p;
		p_ave->update(tn, vn, T_act, T_last, T_sum);
	}


	//saving tissue O2 concentration
	if(do_lum_PlasmaO2_transport&&do_lum_HB_sat_transport&&do_lum_RBC_transport){
        tissueO2_save.push_back(tissueO2s);
    }

    //update for metabolic response
    if(do_metabolic_res){
    	double tn = time.back();
    	Ct_ave->update(tn, tissueO2s, T_act, T_last, T_sum);
    }

    //update CO2 for CO2 control 
    if(do_CO2_control){
    	double tn = time.back();

    	//double PP = CO2_pla_lum->D0_edges[2]->fi[0] ;//co2 concentration
    	double PP = CO2_pla_lum->D0_edges[0]->fi[0] ;//local arterial co2 concentration
    	P_CO2_ave->update(tn, PP, T_act, T_last, T_sum);
    }

}


//--------------------------------------------------------------
void solver_lumped::myogenic_control(double t_act)
{

	// time step
	double dt = t_act - time.back();

	double p = p_ave->average.back();

	// actuator signal
	x_myo = x_myo + dt / tao * (- x_myo + G * (p - p_ref)/(p_ref - atmospheric_pressure/mmHg_to_Pa));
}

//--------------------------------------------------------------
void solver_lumped::set_constants(double g, double rho, double nu, double mmHg, double p0)
{
	gravity = g;
	density = rho;
	kinematic_viscosity = nu;
	mmHg_to_Pa = mmHg;
	atmospheric_pressure = p0;
}

//--------------------------------------------------------------
void solver_lumped::set_non_SI_parameters()
{
	for(int i=0; i<number_of_nodes; i++)
	{
		nodes[i]->pres_ini_non_SI = nodes[i]->pressure_initial/mmHg_to_Pa;
	}
	for(int i=0; i<number_of_edges; i++)
	{
		edges[i]->vfr_ini_non_SI = edges[i]->volume_flow_rate_initial*1.e-6;
		if(edges[i]->type_code == 0) // resistance
		{
			edges[i]->par_non_SI.push_back(edges[i]->parameter[0]/mmHg_to_Pa*1.e-6);
		}
		else if(edges[i]->type_code == 1) // capacitor
		{
			edges[i]->par_non_SI.push_back(edges[i]->parameter[0]*mmHg_to_Pa*1.e6);
		}
		else if(edges[i]->type_code == 2) // elastance
		{
			edges[i]->par_non_SI.push_back(edges[i]->parameter[0]/mmHg_to_Pa*1.e-6);
			edges[i]->par_non_SI.push_back(edges[i]->parameter[1]/mmHg_to_Pa*1.e-6);
		}
		else if(edges[i]->type_code == 3) // inductor
		{
			edges[i]->par_non_SI.push_back(edges[i]->parameter[0]/mmHg_to_Pa*1.e-6);
		}
		else if(edges[i]->type_code == 4) // voltage
		{
			edges[i]->par_non_SI.push_back(edges[i]->parameter[0]/mmHg_to_Pa);
		}
		else if(edges[i]->type_code == 5) // diode
		{
			edges[i]->par_non_SI.push_back(edges[i]->parameter[0]/mmHg_to_Pa*1.e-6);
		}
		else if(edges[i]->type_code == 6) // squared resistance
		{
			edges[i]->par_non_SI.push_back(edges[i]->parameter[0]/mmHg_to_Pa*1.e-6*1.e-6);
		}
		else if(edges[i]->type_code == 7) // resistance_coronary
		{
			edges[i]->par_non_SI.push_back(edges[i]->parameter[0]/mmHg_to_Pa*1.e-6);
		}
		else if(edges[i]->type_code == 8) // capacitor
		{
			edges[i]->par_non_SI.push_back(edges[i]->parameter[0]*mmHg_to_Pa*1.e6);
		}
		else if(edges[i]->type_code == 9) // current source
		{
			edges[i]->par_non_SI.push_back(edges[i]->parameter[0]*1.e6);
		}
		else if(edges[i]->type_code == 10) // piecewise constant resistance
		{
			edges[i]->par_non_SI.push_back(edges[i]->parameter[0]/mmHg_to_Pa*1.e-6);
			edges[i]->par_non_SI.push_back(edges[i]->parameter[1]/mmHg_to_Pa*1.e-6);
		}
	}
}

//--------------------------------------------------------------
void solver_lumped::build_system()
{
	// Clearing the in/out going edges from nodes
	for(unsigned int i=0; i<number_of_nodes; i++)
	{
		nodes[i]->edge_in.clear();
		nodes[i]->edge_out.clear();
	}

	for(unsigned int i=0; i<number_of_edges; i++)
	{
		// starting node
		int node_start = node_id_to_index(edges[i]->node_name_start);

		// ending node
		int node_end = node_id_to_index(edges[i]->node_name_end);

		// saving to edges
		edges[i]->node_index_start = node_start;
		edges[i]->node_index_end = node_end;

		// saving to nodes
		nodes[node_start]->edge_out.push_back(i);
		nodes[node_end]->edge_in.push_back(i);
	}

	set_0D_pointers();

}

//--------------------------------------------------------------
int solver_lumped::node_id_to_index(string node_id)
{
	int i=0, idx=-1;
	bool got_it=false;
	while(i<number_of_nodes && !got_it)
	{
		if(node_id.compare(nodes[i]->name) == 0)
		{
			got_it = true;
			idx = i;
		}
		i++;
	}
	if(idx == -1)
	{
		cout << "\n !!!WARNING!!!\n solver_lumped::node_id_to_index function\nNode is not existing, node_id: " << node_id << "\n Continouing..." << endl;
	}
	return idx;
}

//--------------------------------------------------------------
int solver_lumped::edge_id_to_index(string edge_id)
{
	int i=0, idx=-1;
	bool got_it=false;
	while(i<number_of_edges && !got_it)
	{
		if(edge_id.compare(edges[i]->name) == 0)
		{
			got_it = true;
			idx = i;
		}
		i++;
	}
	if(idx == -1)
	{
		cout << "\n!!!WARNING!!!\n solver_lumped::edge_id_to_index function\nEdge is not existing, edge_id: " << edge_id << "\n Continouing..." << endl;
	}
	return idx;
}

//--------------------------------------------------------------
double solver_lumped::elastance(double t)
{
	vector<double> par{elastance_max_nom,elastance_min_nom};
	return elastance(t,par);
}

//--------------------------------------------------------------
double solver_lumped::elastance(double t, vector<double> par)
{	/*
	// normalized version
	double tn = t * heart_rate/60.;

	// making the elastance periodic
	while(tn>1.)
	{
		tn -= 1.;
	}*/

	//changing heart rate
	double tn = (t - T_sum)*heart_rate/60.;

	double En = 17.4073 * pow(tn,1.9) / (1.+11.2305*pow(tn,1.9)) * 1. / (1.+1.6658e7*pow(tn,21.9));

	//double En = 1.55*pow(tn/0.7,1.9)/(1.+pow(tn/0.7,1.9)) * (1./(1.+pow(tn/1.17,21.9)));

	double E = (par[0]-par[1])*En + par[1];

	// E = E*mmHg_to_Pa*1.e6; // mmHg/ml to SI: Pa/m3 

	return E;
}

//--------------------------------------------------------------
double solver_lumped::elastance_derived(double t, vector<double> par)
{
	/*
	// normalized version
	double tn = t * heart_rate/60.;

	// making the elastance periodic
	while(tn>1.)
	{
		tn -= 1.;
	}*/

	//changing heart rate
	double tn = (t - T_sum)*heart_rate/60.;

	double Enp = (9.450202509727443e-16*pow(tn,0.9) - 1.6570681411267346e-7*pow(tn,22.8) - 2.037762561602155e-6*pow(tn,24.7))/(pow(0.0890432 + pow(tn,1.9),2.)*pow(6.003121623244087e-8 + pow(tn,21.9),2.));

	double Ep = (par[0]-par[1])*Enp;

	//Ep = Ep*mmHg_to_Pa*1.e6;

	return Ep;
}

//--------------------------------------------------------------
void solver_lumped::clear_save_memory()
{
	for(int i=0; i<number_of_nodes; i++)
	{
		nodes[i]->do_save_memory = false;
	}
	for(int i=0; i<number_of_edges; i++)
	{
		edges[i]->do_save_memory = false;
	}
}

//--------------------------------------------------------------
void solver_lumped::set_save_memory(vector<string> edge_list, vector<string> node_list)
{
	for(int i=0; i<edge_list.size(); i++)
	{
		int idx = edge_id_to_index(edge_list[i]);
		if(idx>-1)
		{
			edges[idx]->do_save_memory = true;
		}
	}
	for(int i=0; i<node_list.size(); i++)
	{
		int idx = node_id_to_index(node_list[i]);
		if(idx>-1)
		{
			nodes[idx]->do_save_memory = true;
		}
	}
}


//--------------------------------------------------------------
int NX(double L,double dx, int N) {
    if (floor(L / dx) + 1 > N) {
        return N;
    }
    else {
        return floor(L / dx) + 1;
    }
}


//------------------------------------------------------------
D0_transport::D0_transport(TransportType TType): TType(TType) {
	D0_edges.clear();
}


//------------------------------------------------------------
void D0_transport::update_fi(double dt, solver_lumped& lum_mod, double t_act){
	update_nodes( lum_mod);
	update_edges( dt, lum_mod, t_act);
}


//------------------------------------------------------------
void D0_transport::prescribe_node_fi(TransportType TType, double& finode){
	switch(TType){
	case RBC:
	finode = 4.9e15; // [cell/m3]
	break;

	case C_Plasma_O2:
	finode = 2.9545e-3; // [m3/m3]
	break;

	case HB_O2_saturation:
	finode = 0.97; // [1]
	break;
	} 
}



//--------------------------------------------------------------
void D0_transport::save_variables(){

	for(int i=0;i<D0_edges.size();i++){

		D0_edges[i]->save();
	}

}


//--------------------------------------------------------------
void D0_transport::save_results(string fn, const vector<double>& time, string model_name){
	string file_name, tname;
	switch(this->TType){
	case RBC:
		tname = "RBC";
		break;

	case HB_O2_saturation:
		tname = "HB_O2";
		break;

	case C_Plasma_O2:
		tname = "C_Plasma_O2";
		break;

	case CO2_pla:
		tname = "CO2_pla";
		break;

	case CO2_rbc:
		tname = "CO2_rbc";
		break;

	case HCO3_pla:
		tname = "HCO3_pla";
		break;

	case HCO3_rbc:
		tname = "HCO3_rbc";
		break;

	case HbCO2:
		tname = "HbCO2";
		break;
	}

	mkdir(("results/" + fn + "/" + model_name).c_str(),0777);
	mkdir(("results/" + fn + "/" + model_name + "/" + tname).c_str(),0777);

	for(int i=0;i< D0_edges.size();i++){
		file_name = "results/" + fn + "/" + model_name + "/" + tname + "/" + D0_edges[i]->D0_name + ".txt";
		save_vector(file_name, D0_edges[i]->fi_start, D0_edges[i]->fi_end, time);
	}
	
}

//--------------------------------------------------------------
void D0_transport::save_vector(string folder_name, const vector<double>& st, const vector<double>& en, const vector<double>& time){
    FILE *out_file = fopen(folder_name.c_str(),"w");

	for(unsigned int j=0; j<st.size(); j++)
	{
		double t = time[j];
		double fi_start = st[j];
		double fi_end = en[j];
		fprintf(out_file, "%9.7e, %9.7e, %9.7e\n", t, fi_start, fi_end);
	}
    fclose(out_file);
}

//--------------------------------------------------------------
void D0_transport::save_vector(string folder_name, const vector<double>& vect, const vector<double>& time){
    FILE *out_file = fopen(folder_name.c_str(),"w");

	for(unsigned int j=0; j < time.size(); j++)
	{
		double t = time[j];
		double fi = vect[j];
		fprintf(out_file, "%9.7e, %9.7e\n", t, fi);
	}
    fclose(out_file);
}

//--------------------------------------------------------------
void D0_transport::set_save_memory(){
    do_save_results = true;
    for(int i=0;i<D0_edges.size();i++){
    	D0_edges[i]->do_save_memory=true;
    }
};

//--------------------------------------------------------------
void solver_lumped::O2transport(double dt){

	
	//only one capillary is allowed per 0D model.
	vector<double> HBold = per_cap_BH->fi;
	vector<double>& HB = per_cap_BH->fi;

	vector<double>& plasmaO2 = per_cap_PO2->fi;
	vector<double> plasmaO2old = per_cap_PO2->fi;

	vector<double> RBC = per_cap_RBC->fi;
	vector<double> tissueO2vold = tissueO2v;

	//if(name=="p10"){
	//Mmax = 2.5e-4*2.0;
//}

	//vfr_edge, A, nx are the same for these D0_edges
	int n = per_cap_RBC-> nx;
	double dx = per_cap_RBC-> dx;
	double v = per_cap_RBC->corr_edge->vfr/per_cap_RBC->A*ml_to_m3;

	//BCs
	double fiStartNodePlasma = per_cap_PO2->node_start->PlasmaO2_0Dn;
	double fiStartNodeHB = per_cap_BH->node_start->HBsat_0Dn;
	double fiEndNodePlasma = per_cap_PO2->node_end->PlasmaO2_0Dn;
	double fiEndNodeHB = per_cap_BH->node_end->HBsat_0Dn;

	//capillary plasma concentration
    for (int i = 1; i < n - 1; i++) {
        double Cc1der;
        double HB1der;

        if (v > 0.) {
            Cc1der = (plasmaO2old[i] - plasmaO2old[i-1])/dx;
            HB1der = (HBold[i] - HBold[i-1])/dx;
        }
        else {
            Cc1der = (plasmaO2old[i+1] - plasmaO2old[i])/dx;
            HB1der = (HBold[i+1] - HBold[i])/dx;
        }
        double DCO2 = dCO2_plasma(plasmaO2old[i], HBold[i] , RBC[i]);
        plasmaO2[i] = plasmaO2old[i] + dt*(-v*Cc1der - kc/hc*S_V_c*(plasmaO2old[i]/alpha_b - tissueO2vold[i]/alpha_t) + DCO2/taoO2);
        HB[i] = (1-dt/taoO2)*HBold[i] + dt/taoO2*HBsat_equilibrium(plasmaO2old[i] / alpha_b) - dt*v*HB1der;
    }

    // plasma and HBsat BC
    if (v > 0.) {
        double Cc1der = (plasmaO2old[n - 1] - plasmaO2old[n - 2])/dx;
        double HB1der = (HBold[n - 1] - HBold[n - 2])/dx;
        double DCO2 = dCO2_plasma(plasmaO2old[n - 1], HBold[n - 1] , RBC[n - 1]);
        plasmaO2[n - 1] = plasmaO2old[n - 1] + dt*(-v*Cc1der - kc/hc*S_V_c*(plasmaO2old[n-1]/alpha_b - tissueO2vold[n-1]/alpha_t) + DCO2 / taoO2);
        plasmaO2[0] = fiStartNodePlasma;
        HB[n-1] = (1-dt/taoO2)*HBold[n-1] + dt/taoO2*HBsat_equilibrium(plasmaO2old[n-1] / alpha_b) - dt*v*HB1der;
        HB[0] = fiStartNodeHB;
    }
    else {
    	double Cc1der = (plasmaO2old[1] - plasmaO2old[0])/dx;
    	double HB1der = (HBold[1] - HBold[0])/dx;
    	double DCO2 = dCO2_plasma(plasmaO2old[0], HBold[0] , RBC[0]);
        plasmaO2[n - 1] = fiEndNodePlasma; 
        plasmaO2[0] = plasmaO2old[0] + dt*(-v*Cc1der - kc/hc*S_V_c*(plasmaO2old[0]/alpha_b - tissueO2vold[0]/alpha_t) + DCO2/taoO2);
        HB[n-1] = fiEndNodeHB;
        HB[0] = (1-dt/taoO2)*HBold[0] + dt/taoO2*HBsat_equilibrium(plasmaO2old[0] / alpha_b)- dt*v*HB1der;
    }

    //tissue concentration
    for(int i=0; i<n; i++){
    tissueO2v[i] = tissueO2vold[i] + dt*( vessel_dilation(1) * kc/hc*S_V_c*fi_c/fi_t*(plasmaO2old[i]/alpha_b - tissueO2vold[i]/alpha_t) - Mmax*tissueO2vold[i]/(tissueO2vold[i] + C50 ));

    }

    tissueO2s = average(tissueO2v);
};


//--------------------------------------------------------------
void solver_lumped::pulmonary_O2transport(double dt){
	double v = pul_cap_BH->corr_edge->vfr/pul_cap_BH->A * ml_to_m3;
	//cout<<v<<endl;

	vector<double> HBold = pul_cap_BH-> fi;
	vector<double>& HB = pul_cap_BH-> fi;

	vector<double>& plasmaO2 = pul_cap_PO2-> fi;
	vector<double> plasmaO2old = pul_cap_PO2-> fi;

	vector<double> RBC = pul_cap_RBC-> fi;

	int n = pul_cap_RBC-> nx;
	vector<double> K_pul_v = sin_2(K_pul_scale, n);

	//vfr_edge, A, nx are the same for these D0_edges
	double dx = pul_cap_RBC-> dx;

	//BCs
	double fiStartNodePlasma = pul_cap_PO2->node_start->PlasmaO2_0Dn;
	double fiStartNodeHB = pul_cap_BH->node_start->HBsat_0Dn;
	double fiEndNodePlasma = pul_cap_PO2->node_end->PlasmaO2_0Dn;
	double fiEndNodeHB = pul_cap_BH->node_end->HBsat_0Dn;

    for (int i = 1; i < n - 1; i++) {
        double Cc1der;
        double HB1der;

        if (v > 0.) {
            Cc1der = (plasmaO2old[i] - plasmaO2old[i-1])/dx;
            HB1der = (HBold[i] - HBold[i-1])/dx;
        }
        else {
            Cc1der = (plasmaO2old[i+1] - plasmaO2old[i])/dx;
            HB1der = (HBold[i+1] - HBold[i])/dx;
        }
        double DCO2 = dCO2_plasma(plasmaO2old[i], HBold[i] , RBC[i]);
        plasmaO2[i] = plasmaO2old[i] + dt*(-v*Cc1der - K_pul_v[i] *(plasmaO2old[i]/alpha_b - PO2_alveolar) + DCO2/taoO2_p);
        HB[i] = (1-dt/taoO2_p)*HBold[i] + dt/taoO2_p*HBsat_equilibrium(plasmaO2old[i] / alpha_b) - dt*v*HB1der;
    }

    // plasma and HBsat BC
    if (v > 0.) {
        double Cc1der = (plasmaO2old[n - 1] - plasmaO2old[n - 2])/dx;
        double HB1der = (HBold[n - 1] - HBold[n - 2])/dx;
        double DCO2 = dCO2_plasma(plasmaO2old[n - 1], HBold[n - 1] , RBC[n - 1]);
        plasmaO2[n - 1] = plasmaO2old[n - 1] + dt*(-v*Cc1der - K_pul_v[n-1] *(plasmaO2old[n-1]/alpha_b - PO2_alveolar) + DCO2 / taoO2_p);
        plasmaO2[0] = fiStartNodePlasma;
        HB[n-1] = (1-dt/taoO2_p)*HBold[n-1] + dt/taoO2_p*HBsat_equilibrium(plasmaO2old[n-1] / alpha_b) - dt*v*HB1der;
        HB[0] = fiStartNodeHB;
    }
    else {
    	double Cc1der = (plasmaO2old[1] - plasmaO2old[0])/dx;
    	double HB1der = (HBold[1] - HBold[0])/dx;
    	double DCO2 = dCO2_plasma(plasmaO2old[0], HBold[0] , RBC[0]);
        plasmaO2[n - 1] = fiEndNodePlasma; 
        plasmaO2[0] = plasmaO2old[0] + dt*(-v*Cc1der - K_pul_v[0] *(plasmaO2old[0]/alpha_b - PO2_alveolar) + DCO2/taoO2_p);
        HB[n-1] = fiEndNodeHB;
        HB[0] = (1-dt/taoO2_p)*HBold[0] + dt/taoO2_p*HBsat_equilibrium(plasmaO2old[0] / alpha_b)- dt*v*HB1der;
    }
    //cout<<HB[0]<<endl;

    //for(int i=0; i<n; i++){
    //	cout<<HB[i]<<" ";
    //}
    //cout<<endl<<endl;


};


//--------------------------------------------------------------
double solver_lumped::dCO2_plasma(double CO2_plasma_old, double HBsat_old, double C_RBC){
    double HBsat_eq = HBsat_equilibrium(CO2_plasma_old/alpha_b);
    double dHBsat = HBsat_old - HBsat_eq;
    return dHBsat*C_RBC*Z;
}


//--------------------------------------------------------------
double solver_lumped::HBsat_equilibrium(double PO2){
	return L_HBsat/(1. + exp(-k_HBsat*( PO2 - m_HBsat ))) + b_HBsat;
}


//--------------------------------------------------------------
void solver_lumped::init_lum_tissueO2(){
    tissueO2v.clear();
    tissueO2_save.clear();

    //only one peripheral capillary is allowed in a lumped model
    int t;

	for(int i=0;i<HBsatlum ->D0_edges.size();i++){
		if(HBsatlum ->D0_edges[i]->is_per_capillary){
			t=HBsatlum ->D0_edges[i]->nx;
    		//O2 transport initialization
    		tissueO2s = init_tissueO2;
    		tissueO2v.assign( t , init_tissueO2);
		}
	}
}


//--------------------------------------------------------------
double solver_lumped::turn_source(double t){
	return 1/(1 + exp(-(t-10.)));
}


//--------------------------------------------------------------
void solver_lumped::save_tissueO2(string folder_name, const vector<double>& time){
		string file_name = folder_name;
		HBsatlum->save_vector(file_name, tissueO2_save, time);
}


//--------------------------------------------------------------
vector<double> D0_transport::linear_dist(double avg, double dist, int len){
	vector<double> r;
	r.assign(len, 0.);
	double d = dist / (len - 1);
	int z = len-1;
	if (len%2 == 1) {
		for (int i = -((len - 1) / 2); i < (len - 1) / 2 + 1; i++) {
			r[z] = avg + i * d;
			z--;
		}
	}
	else {
		for (int i = -(len / 2); i < len / 2 ; i++) {
			r[z] = avg + (i + 0.5) * d;
			z--;
		}
	}
	return r;
}


//--------------------------------------------------------------
vector<double> solver_lumped::sin_2(double scale, int nx){
	vector<double> r;
	r.assign(nx, 0.);
	for(int i=0; i<nx; i++){
		r[i] = sin( i*pi/(nx-1) )*sin( i*pi/(nx-1) )*scale;
	}
	return r;
}


//--------------------------------------------------------------
void solver_lumped::assign_perif_O2_params(vector<string> sv){
    fi_c = stod(sv[1],0);
    fi_t = stod(sv[2],0);
    alpha_b = stod(sv[3],0);
    alpha_t = stod(sv[4],0);
    hc = stod(sv[5],0);
    S_V_c = stod(sv[6],0);
    kc = stod(sv[7],0);
    //Mmax = stod(sv[8],0); // Mmax is read from the model files 
    C50 = stod(sv[9],0);
    taoO2 = stod(sv[10],0);
    Z = stod(sv[11],0);
}


//--------------------------------------------------------------
void solver_lumped::assign_haemogobin_sat_params(vector<string> sv){
	L_HBsat = stod(sv[1],0);
	k_HBsat = stod(sv[2],0);
	b_HBsat = stod(sv[3],0);
	m_HBsat = stod(sv[4],0);
}


//--------------------------------------------------------------
void solver_lumped::assign_pulmonary_O2_params(vector<string> sv){
	PO2_alveolar = stod(sv[1],0);
	//K_pul_O2 = stod(sv[2],0);
	taoO2_p = stod(sv[3],0);
	K_pul_scale =stod(sv[4],0) ;
}


//--------------------------------------------------------------
void solver_lumped::metabolic_response(double t_act)
{

	// time step
	double dt = t_act - time.back();

	double Ct = Ct_ave->average.back();

	// actuator signal
	x_met = x_met + dt / tao_met * (- x_met + G_met * (Ct - Ct_ref)/Ct_ref);

}


double solver_lumped::vessel_dilation(int edgeindex){
	//return 1.;
	return pow( edges[edgeindex]->parameter_factor ,-0.5);
}


//--------------------------------------------------------------------------------------------------
D0_edge::D0_edge(string D0_name, double L, double A, int  nx, TransportType TType, double init, string node_s_name, string node_e_name, string corr_edge_name):D0_name(D0_name),L(L),A(A),nx(nx),TType(TType),init(init),node_s_name(node_s_name),node_e_name(node_e_name),corr_edge_name(corr_edge_name){

	dx = L/(nx-1);
	fi_start.clear();
	fi_end.clear();
	fi.clear();

	fi.assign(nx, init);
	if(do_save_memory){
		save();}

}


//--------------------------------------------------------------------------------------------------
void D0_transport::update_edges( double dt, solver_lumped& lum_mod, double t_act){
//capillary edges are updates from solver_lumped
	for(int i=0; i< D0_edges.size(); i++ ){

		if((!D0_edges[i]->is_per_capillary && !D0_edges[i]->is_pul_capillary)){
			if(D0_edges[i]->is_diode){
				D0_edges[i]->update_diode();
			}
			else if(D0_edges[i]->is_capacitor){
				D0_edges[i]->update_capacitor(dt);
			}
			else if(D0_edges[i]->is_elastance){
				double E = lum_mod.elastance(t_act, D0_edges[i]->corr_edge->parameter);
				E = E*mmHg_to_Pa*1.e6; // mmHg/ml to SI: Pa/m3
				D0_edges[i]->update_elastance(dt, E);
			}
			else{
				D0_edges[i]->virt1D(dt);
			}
		}
	}
}


//--------------------------------------------------------------------------------------------------
void D0_edge::save(){

	fi_start.push_back(fi[0]);
	fi_end.push_back(fi.back());
}


//--------------------------------------------------------------------------------------------------
void solver_lumped::set_0D_pointers(){

	if(do_lum_RBC_transport){
		for(int i=0; i<RBClum->D0_edges.size(); i++){
			for(int j=0; j<edges.size();j++){
				if(RBClum->D0_edges[i]->corr_edge_name == edges[j]->name ){RBClum->D0_edges[i]->corr_edge = edges[j];}
			}
		}
	}


	if(do_lum_PlasmaO2_transport){
		for(int i=0; i<PlasmaO2lum->D0_edges.size(); i++){
			for(int j=0; j<edges.size();j++){
				if(PlasmaO2lum->D0_edges[i]->corr_edge_name == edges[j]->name ){PlasmaO2lum->D0_edges[i]->corr_edge = edges[j];}
			}
		}
	}

	if(do_lum_HB_sat_transport){
		for(int i=0; i<HBsatlum->D0_edges.size(); i++){
			for(int j=0; j<edges.size();j++){
				if(HBsatlum->D0_edges[i]->corr_edge_name == edges[j]->name ){HBsatlum->D0_edges[i]->corr_edge = edges[j];}
			}
		}
	}

	//CO2 transport
	if(do_lum_pla_CO2_transport){
		for(int i=0; i<CO2_pla_lum->D0_edges.size(); i++){
			for(int j=0; j<edges.size();j++){
				if(CO2_pla_lum->D0_edges[i]->corr_edge_name == edges[j]->name ){CO2_pla_lum->D0_edges[i]->corr_edge = edges[j];}
			}
		}
	}

	if(do_lum_rbc_CO2_transport){
		for(int i=0; i<CO2_rbc_lum->D0_edges.size(); i++){
			for(int j=0; j<edges.size();j++){
				if(CO2_rbc_lum->D0_edges[i]->corr_edge_name == edges[j]->name ){CO2_rbc_lum->D0_edges[i]->corr_edge = edges[j];}
			}
		}
	}

	if(do_lum_pla_HCO3_transport){
		for(int i=0; i<HCO3_pla_lum->D0_edges.size(); i++){
			for(int j=0; j<edges.size();j++){
				if(HCO3_pla_lum->D0_edges[i]->corr_edge_name == edges[j]->name ){HCO3_pla_lum->D0_edges[i]->corr_edge = edges[j];}
			}
		}
	}

	if(do_lum_rbc_HCO3_transport){
		for(int i=0; i<HCO3_rbc_lum->D0_edges.size(); i++){
			for(int j=0; j<edges.size();j++){
				if(HCO3_rbc_lum->D0_edges[i]->corr_edge_name == edges[j]->name ){HCO3_rbc_lum->D0_edges[i]->corr_edge = edges[j];}
			}
		}
	}

	if(do_lum_HbCO2_transport){
		for(int i=0; i<HbCO2_lum->D0_edges.size(); i++){
			for(int j=0; j<edges.size();j++){
				if(HbCO2_lum->D0_edges[i]->corr_edge_name == edges[j]->name ){HbCO2_lum->D0_edges[i]->corr_edge = edges[j];}
			}
		}
	}


}


//--------------------------------------------------------------------------------------------------
void solver_lumped::capillary_O2_transport(double dt){

	if(do_pul_O2_rtansport){

		pulmonary_O2transport( dt);
	}

    if(do_per_O2_rtansport){

    	O2transport( dt);
    }

};


//--------------------------------------------------------------------------------------------------
void solver_lumped::autoregulation(double t_act){

	if(t_act<=80.){return;}

	if(do_myogenic)
	{
		//updates x_myo
		myogenic_control(t_act);
	}

	if(do_metabolic_res){
		//updates x_met
		metabolic_response(t_act);
	}

	if(do_CO2_control){
		//updates x_met
		CO2_response(t_act);
	}

	if(do_perif_baroreflex){
		//updates x_bar
		perif_baroreflex();
	}



	//updates the parameter factor of the resistance
	if(do_CO2_control || do_metabolic_res || do_myogenic || do_perif_baroreflex){
	update_R_fact();}
}

/*
//--------------------------------------------------------------------------------------------------
void solver_lumped::update_R_fact(){

	vector<int> Ridx{0,1}; // which resistors are we modifying

	double FF; 
	for(int i=0; i<Ridx.size(); i++)
	{
		double Rmax, Rmin; // calculated from r_min, r_max from "Regulation of Coronary Microvascular Resistance in Health and Disease" pic 12.2
		if((x_met + x_myo - x_CO2) < 0){ // x_co2 has the opposite effect -> negative sign
			//the sigmoid curve is the same for the two responses
			Rmin = sat1_met;
			Rmax = 2. - Rmin;
		}
		else{
			Rmax = sat2_met;
			Rmin = 2. - sat2_met;
		}

		double ff = 80. / (Rmax - Rmin) ;
		FF = (Rmax + Rmin * exp(- (x_met + x_myo - x_CO2) * ff)) / (1. + exp(-(x_met + x_myo - x_CO2) * ff));
		

		if(do_CO2_control){
			//cout<<x_met<<"  "<<x_myo<<"  "<<x_CO2<<endl;
			cout.precision(10);
			cout << FF << endl<<endl;}

		edges[Ridx[i]]->parameter_factor = FF;
	}

}*/


//--------------------------------------------------------------------------------------------------
void solver_lumped::update_R_fact(){

	vector<int> Ridx{0,1}; // which resistors are we modifying


	double Rmax = sat2_met;
	double Rmin = sat1_met;
	double x = x_met + x_myo - x_CO2 + x_bar;
	double A = (Rmax - Rmin) / (1.0 - Rmin) - 1.0;
	double FF = Rmin + (Rmax - Rmin) / (1.0 + A * exp(- 30. * x));


	//cout.precision(10);
	cout << FF << endl;

	for(int i=0; i<Ridx.size(); i++)
	{
		edges[Ridx[i]]->parameter_factor = FF;
	}

	
 


}




//--------------------------------------------------------------------------------------------------
void D0_edge::update_capacitor(double dt){
	double V = corr_edge->parameter[0] * abs( node_start->p - node_end->p  )/mmHg_to_Pa; // SI
	double fi_old = fi[0];
	//cout<<( node_start->p - node_end->p  )<<endl;


	double Q = corr_edge->vfr;

	//if( (D0_name == "c_pa") && (Q>0.) && (node_start->PlasmaO2_0Dn != 0.0011467116) ){cout<< "Alma";}
	//if(D0_name == "c_pa"){cout<<Q<<endl<<node_start->PlasmaO2_0Dn<<endl<<endl;}


	if(Q>0.){
		double f;
		switch(TType){
			case RBC:
				f=node_start->RBC_fi0Dn;
			break;

			case HB_O2_saturation:
				f=node_start->HBsat_0Dn;
			break;

			case C_Plasma_O2:
				f=node_start->PlasmaO2_0Dn;
			break;

			//CO2
			case CO2_pla:
				f=node_start->CO2_pla_n;
			break;

			case CO2_rbc:
				f=node_start->CO2_rbc_n;
			break;

			case HCO3_pla:
				f=node_start->HCO3_pla_n;
			break;

			case HCO3_rbc:
				f=node_start->HCO3_rbc_n;
			break;

			case HbCO2:
				f=node_start->HbCO2_n;
			break;
		}
		if(V+Q*dt*ml_to_m3 != 0.){
		fi[0] = (fi_old*V + Q*dt*f*ml_to_m3)/(V+Q*dt*ml_to_m3);}

		//if ((D0_name == "c_pa") && (Q > 0.) && (node_start->PlasmaO2_0Dn != 0.0011467116)) {cout << std::setprecision(17) << V << endl << Q*dt*ml_to_m3 <<endl<<endl;}
		//if ((D0_name == "c_pa") && (Q > 0.) && (node_start->PlasmaO2_0Dn != 0.0011467116)) {cout << std::setprecision(17) << fi[0] << endl << endl;}

	}

}




//--------------------------------------------------------------------------------------------------
void D0_edge::update_elastance(double dt, double E){
	double Q = corr_edge->vfr;
	double dV = Q*dt*ml_to_m3;
	double V = V0 + abs( node_end->p - node_start->p  )/E/mmHg_to_Pa + dV; //SI

	double fi_old = fi[0];



	if(Q<0.){
		double f;
		switch(TType){
			case RBC:
				f=node_end->RBC_fi0Dn;
			break;

			case HB_O2_saturation:
				f=node_end->HBsat_0Dn;
			break;

			case C_Plasma_O2:
				f=node_end->PlasmaO2_0Dn;
			break;

			//CO2
			case CO2_pla:
				f=node_end->CO2_pla_n;
			break;

			case CO2_rbc:
				f=node_end->CO2_rbc_n;
			break;

			case HCO3_pla:
				f=node_end->HCO3_pla_n;
			break;

			case HCO3_rbc:
				f=node_end->HCO3_rbc_n;
			break;

			case HbCO2:
				f=node_end->HbCO2_n;
			break;

		}

		if((V-dV)!=0){
		fi[0] = (fi_old*V - dV*f)/(V-dV);}
		//cout<<f<<"  "<<fi_old<<"  "<<fi[0]<<endl;
	}

}

//--------------------------------------------------------------------------------------------------
bool solver_lumped::check_whole_period(double t_act){
	//sum of the previous periods
	if(T_sum + T_act < t_act){ //passes the end of the next period in this instant
		return true;
	}

	return false;
}

//--------------------------------------------------------------------------------------------------
void solver_lumped::update_period_time(double T_act_new){
	T_sum += T_act; // end of the last completed cycle
	T_last = T_act;
	T_act = T_act_new;
	heart_rate = 60./T_act_new;
	time_period = T_act_new;

	//new period
	period ++;
}

/*
----------------------------------------------------------------------------------------------------
CO2 transport
----------------------------------------------------------------------------------------------------
*/

//--------------------------------------------------------------------------------------------------
double solver_lumped::eta_hco3(double HCO3_rbc, double HCO3_pla, double CO2_pla, double CO2_rbc){
	double g3 = CO2_pla*fi_pla + CO2_rbc*fi_rbc;
	double r_hco3 = (462*exp(6.2255*g3)-340*exp(-66.7556*g3)+0.62)*0.95*1e-3 - g3;
	//cout<<r_hco3<<endl;
	//cout<<g3<<endl;
	return r_hco3 - (HCO3_rbc*fi_rbc+HCO3_pla*fi_pla);
};

//--------------------------------------------------------------------------------------------------
double solver_lumped::eta_hb(double HbCO2, double CO2_pla, double CO2_rbc){
	double g3 = CO2_pla*fi_pla + CO2_rbc*fi_rbc;
	//cout<<g3<<endl;
	double r_hb = (462*exp(6.2255*g3*fi_pla)-340*exp(-66.7556*g3*fi_pla)+0.62)*0.05*1e-3;
 	return r_hb - HbCO2*fi_rbc;
};

//--------------------------------------------------------------------------------------------------
void solver_lumped::CO2transport(double dt){
	//only one capillary is allowed per 0D model.
	vector<double> CO2_pla_old = per_cap_CO2_pla->fi;
	vector<double>& CO2_pla = per_cap_CO2_pla->fi;

	vector<double> CO2_rbc_old = per_cap_CO2_rbc->fi;
	vector<double>& CO2_rbc = per_cap_CO2_rbc->fi;

	vector<double> HCO3_pla_old = per_cap_HCO3_pla->fi;
	vector<double>& HCO3_pla = per_cap_HCO3_pla->fi;

	vector<double> HCO3_rbc_old = per_cap_HCO3_rbc->fi;
	vector<double>& HCO3_rbc = per_cap_HCO3_rbc->fi;

	vector<double> HbCO2_old = per_cap_HbCO2->fi;
	vector<double>& HbCO2 = per_cap_HbCO2->fi;

	vector<double> tissueCO2vold = tissueCO2v;

	//vfr_edge, A, nx are the same for these D0_edges
	int n = per_cap_CO2_pla-> nx;
	double dx = per_cap_CO2_pla-> dx;
	double v = per_cap_CO2_pla->corr_edge->vfr/per_cap_CO2_pla->A*ml_to_m3;

	//BCs
	double CO2_pla_n_s = per_cap_CO2_pla->node_start->CO2_pla_n; //node start
	double CO2_pla_n_e = per_cap_CO2_pla->node_end->CO2_pla_n; //node end
	double CO2_rbc_n_s = per_cap_CO2_rbc->node_start->CO2_rbc_n;
	double CO2_rbc_n_e = per_cap_CO2_rbc->node_end->CO2_rbc_n;
	double HCO3_pla_n_s = per_cap_HCO3_pla->node_start->HCO3_pla_n;
	double HCO3_pla_n_e = per_cap_HCO3_pla->node_end->HCO3_pla_n;
	double HCO3_rbc_n_s = per_cap_HCO3_rbc->node_start->HCO3_rbc_n;
	double HCO3_rbc_n_e = per_cap_HCO3_rbc->node_end->HCO3_rbc_n;
	double HbCO2_n_s = per_cap_HbCO2->node_start->HbCO2_n;
	double HbCO2_n_e = per_cap_HbCO2->node_end->HbCO2_n;


	//capillary
    for (int i = 1; i < n - 1; i++) {
        double CO2_pla_der;
        double CO2_rbc_der;
        double HCO3_pla_der;
        double HCO3_rbc_der;
        double HbCO2_der;

        if (v > 0.) {
            CO2_pla_der = (CO2_pla_old[i] - CO2_pla_old[i-1])/dx;
            CO2_rbc_der = (CO2_rbc_old[i] - CO2_rbc_old[i-1])/dx;
            HCO3_pla_der = (HCO3_pla_old[i] - HCO3_pla_old[i-1])/dx;
            HCO3_rbc_der = (HCO3_rbc_old[i] - HCO3_rbc_old[i-1])/dx;
            HbCO2_der = (HbCO2_old[i] - HbCO2_old[i-1])/dx;
        }
        else {
            CO2_pla_der = (CO2_pla_old[i+1] - CO2_pla_old[i])/dx;
            CO2_rbc_der = (CO2_rbc_old[i+1] - CO2_rbc_old[i])/dx;
            HCO3_pla_der = (HCO3_pla_old[i+1] - HCO3_pla_old[i])/dx;
            HCO3_rbc_der = (HCO3_rbc_old[i+1] - HCO3_rbc_old[i])/dx;
            HbCO2_der = (HbCO2_old[i+1] - HbCO2_old[i])/dx;
        }


        double Eta_hco3 = eta_hco3(HCO3_rbc_old[i], HCO3_pla_old[i], CO2_pla_old[i], CO2_rbc_old[i]);
        double Eta_hb = eta_hb( HbCO2_old[i], CO2_pla_old[i], CO2_rbc_old[i]);


        //CO2_pla[i] = CO2_pla_old[i] - dt/(fi_pla*fi_c)*v*CO2_pla_der*ksi_c*ksi_pla - dt/tao_co2_pla_rbc*(CO2_pla_old[i] - CO2_rbc_old[i]*alpha_co2_pla/alpha_co2_rbc) + dt/tao_co2_pla_tis*(tissueCO2vold[i]*alpha_co2_pla/alpha_co2_tis - CO2_pla_old[i]);
        CO2_pla[i] = CO2_pla_old[i] - dt/(fi_pla*fi_c)*v*CO2_pla_der*ksi_c*ksi_pla - dt/tao_co2_pla_rbc*(CO2_pla_old[i] - CO2_rbc_old[i]*alpha_co2_pla/alpha_co2_rbc) + dt*Mmax*tissueO2v[i]/(tissueO2v[i]+C50)*RQ * fi_t/fi_c/fi_pla;

        CO2_rbc[i] = CO2_rbc_old[i] - dt/fi_rbc*v*CO2_rbc_der*ksi_rbc + dt/tao_co2_rbc_pla*(CO2_pla_old[i]*alpha_co2_rbc/alpha_co2_pla - CO2_rbc_old[i]) - dt/tao_hco3*Eta_hco3 - dt/tao_hbco2*Eta_hb;
        HCO3_pla[i] = HCO3_pla_old[i] - dt/fi_pla*v*HCO3_pla_der*ksi_pla - dt/tao_hco3_pla_rbc*(HCO3_pla_old[i]-HCO3_rbc_old[i]/alpha_hco3_rbc*alpha_hco3_pla);
        HCO3_rbc[i] = HCO3_rbc_old[i] - dt/fi_rbc*v*HCO3_rbc_der*ksi_rbc + dt/tao_hco3_rbc_pla*(HCO3_pla_old[i]*alpha_hco3_rbc/alpha_hco3_pla-HCO3_rbc_old[i]) + dt/tao_hco3*Eta_hco3;
        HbCO2[i] = HbCO2_old[i] - v*dt*HbCO2_der + dt/tao_hbco2*Eta_hb;
    }

    //BCs
    if (v > 0.) {
        double CO2_pla_der = (CO2_pla_old[n - 1] - CO2_pla_old[n - 2])/dx;
        double CO2_rbc_der = (CO2_rbc_old[n - 1] - CO2_rbc_old[n - 2])/dx;
        double HCO3_pla_der = (HCO3_pla_old[n - 1] - HCO3_pla_old[n - 2])/dx;
        double HCO3_rbc_der = (HCO3_rbc_old[n - 1] - HCO3_rbc_old[n - 2])/dx;
        double HbCO2_der = (HbCO2_old[n - 1] - HbCO2_old[n - 2])/dx;

        double Eta_hco3 = eta_hco3(HCO3_rbc_old[n - 1], HCO3_pla_old[n - 1], CO2_pla_old[n - 1], CO2_rbc_old[n - 1]);
        double Eta_hb = eta_hb( HbCO2_old[n - 1], CO2_pla_old[n - 1] , CO2_rbc_old[n - 1]);


        //CO2_pla[n-1] = CO2_pla_old[n-1] - dt/(fi_pla*fi_c)*v*CO2_pla_der*ksi_c*ksi_pla - dt/tao_co2_pla_rbc*(CO2_pla_old[n-1] - CO2_rbc_old[n-1]*alpha_co2_pla/alpha_co2_rbc) + dt/tao_co2_pla_tis*(tissueCO2vold[n-1]*alpha_co2_pla/alpha_co2_tis - CO2_pla_old[n-1]);
        CO2_pla[n-1] = CO2_pla_old[n-1] - dt/(fi_pla*fi_c)*v*CO2_pla_der*ksi_c*ksi_pla - dt/tao_co2_pla_rbc*(CO2_pla_old[n-1] - CO2_rbc_old[n-1]*alpha_co2_pla/alpha_co2_rbc) + dt*Mmax*tissueO2v[n-1]/(tissueO2v[n-1]+C50)*RQ* fi_t/fi_c/fi_pla;

        CO2_rbc[n-1] = CO2_rbc_old[n-1] - dt/fi_rbc*v*CO2_rbc_der*ksi_rbc + dt/tao_co2_rbc_pla*(CO2_pla_old[n-1]*alpha_co2_rbc/alpha_co2_pla - CO2_rbc_old[n-1]) - dt/tao_hco3*Eta_hco3 - dt/tao_hbco2*Eta_hb;
        HCO3_pla[n-1] = HCO3_pla_old[n-1] - dt/fi_pla*v*HCO3_pla_der*ksi_pla - dt/tao_hco3_pla_rbc*(HCO3_pla_old[n-1]-HCO3_rbc_old[n-1]/alpha_hco3_rbc*alpha_hco3_pla);
        HCO3_rbc[n-1] = HCO3_rbc_old[n-1] - dt/fi_rbc*v*HCO3_rbc_der*ksi_rbc + dt/tao_hco3_rbc_pla*(HCO3_pla_old[n-1]*alpha_hco3_rbc/alpha_hco3_pla-HCO3_rbc_old[n-1]) + dt/tao_hco3*Eta_hco3;
        HbCO2[n-1] = HbCO2_old[n-1] - v*dt*HbCO2_der+dt/tao_hbco2*Eta_hb;

        CO2_pla[0] = CO2_pla_n_s;
        CO2_rbc[0] = CO2_rbc_n_s;
        HCO3_pla[0] = HCO3_pla_n_s;
        HCO3_rbc[0] = HCO3_rbc_n_s;
        HbCO2[0] = HbCO2_n_s;


    }
    else {
        double CO2_pla_der = (CO2_pla_old[1] - CO2_pla_old[0])/dx;
        double CO2_rbc_der = (CO2_rbc_old[1] - CO2_rbc_old[0])/dx;
        double HCO3_pla_der = (HCO3_pla_old[1] - HCO3_pla_old[0])/dx;
        double HCO3_rbc_der = (HCO3_rbc_old[1] - HCO3_rbc_old[0])/dx;
        double HbCO2_der = (HbCO2_old[1] - HbCO2_old[0])/dx;

        double Eta_hco3 = eta_hco3(HCO3_rbc_old[0], HCO3_pla_old[0], CO2_pla_old[0], CO2_rbc_old[0]);
        double Eta_hb = eta_hb( HbCO2_old[0], CO2_pla_old[0], CO2_rbc_old[0]);
        

        //CO2_pla[0] = CO2_pla_old[0] - dt/(fi_pla*fi_c)*v*CO2_pla_der*ksi_c*ksi_pla - dt/tao_co2_pla_rbc*(CO2_pla_old[0] - CO2_rbc_old[0]*alpha_co2_pla/alpha_co2_rbc) + dt/tao_co2_pla_tis*(tissueCO2vold[0]*alpha_co2_pla/alpha_co2_tis - CO2_pla_old[0]);
        CO2_pla[0] = CO2_pla_old[0] - dt/(fi_pla*fi_c)*v*CO2_pla_der*ksi_c*ksi_pla - dt/tao_co2_pla_rbc*(CO2_pla_old[0] - CO2_rbc_old[0]*alpha_co2_pla/alpha_co2_rbc) + dt*Mmax*tissueO2v[0]/(tissueO2v[0]+C50)*RQ * fi_t/fi_c/fi_pla;//
        
        CO2_rbc[0] = CO2_rbc_old[0] - dt/fi_rbc*v*CO2_rbc_der*ksi_rbc + dt/tao_co2_rbc_pla*(CO2_pla_old[0]*alpha_co2_rbc/alpha_co2_pla - CO2_rbc_old[0]) - dt/tao_hco3*Eta_hco3 - dt/tao_hbco2*Eta_hb;
        HCO3_pla[0] = HCO3_pla_old[0] - dt/fi_pla*v*HCO3_pla_der*ksi_pla - dt/tao_hco3_pla_rbc*(HCO3_pla_old[0]-HCO3_rbc_old[0]/alpha_hco3_rbc*alpha_hco3_pla);
        HCO3_rbc[0] = HCO3_rbc_old[0] - dt/fi_rbc*v*HCO3_rbc_der*ksi_rbc + dt/tao_hco3_rbc_pla*(HCO3_pla_old[0]*alpha_hco3_rbc/alpha_hco3_pla-HCO3_rbc_old[0]) + dt/tao_hco3*Eta_hco3;
        HbCO2[0] = HbCO2_old[0] - v*dt*HbCO2_der+dt/tao_hbco2*Eta_hb;

        CO2_pla[n-1] = CO2_pla_n_e;
        CO2_rbc[n-1] = CO2_rbc_n_e;
        HCO3_pla[n-1] = HCO3_pla_n_e;
        HCO3_rbc[n-1] = HCO3_rbc_n_e;
        HbCO2[n-1] = HbCO2_n_e;

    }

/*
    //tissue concentration
    for(int i=0; i<n; i++){
	tissueCO2v[i] =  tissueCO2vold[i] + dt*Mmax*tissueO2v[i]/(tissueO2v[i]+C50)*RQ - dt/tao_co2_tis_pla*(tissueCO2vold[i] - CO2_pla_old[i]*alpha_co2_tis/alpha_co2_pla);
    }
    tissueCO2s = average(tissueCO2v);

    //cout<<tissueCO2s/alpha_co2_tis<<endl;
    cout<<Mmax*tissueO2v[30]/(tissueO2v[30]+C50)*RQ<<"  "<<1/tao_co2_tis_pla*(tissueCO2vold[30] - CO2_pla_old[30]*alpha_co2_tis/alpha_co2_pla)<< "  " <<name<<endl;
    //cout<<Mmax*tissueO2v[30]/(tissueO2v[30]+C50)<<endl;
    //cout<<tissueO2s<<endl;
    cout<<tissueCO2s/alpha_co2_tis/mmHg_to_Pa<<endl<<endl;*/

    //for (int i = 1; i < n - 1; i++){cout<<CO2_pla[i]<<"  ";}
    //cout<<endl;
    //cout<<CO2_pla[n-1]<<endl;

};


//--------------------------------------------------------------
void solver_lumped::assign_perif_CO2_params(vector<string> sv){
    tao_hco3_rbc_pla = stod(sv[1],0);
 	tao_hco3_pla_rbc = stod(sv[2],0);
 	tao_co2_pla_tis = stod(sv[3],0);
 	tao_co2_tis_pla = stod(sv[4],0);
 	tao_co2_rbc_pla = stod(sv[5],0);
 	tao_co2_pla_rbc = stod(sv[6],0);

 	tao_hco3 = stod(sv[7],0);
 	tao_hbco2 = stod(sv[8],0);

 	RQ = stod(sv[9],0);

 	alpha_co2_rbc = stod(sv[10],0);
    alpha_co2_pla = stod(sv[11],0);
    alpha_co2_tis = stod(sv[12],0);
    alpha_hco3_rbc = stod(sv[13],0);
    alpha_hco3_pla = stod(sv[14],0);

 	ksi_pla = stod(sv[15],0);
 	ksi_rbc = stod(sv[16],0);
 	ksi_c = stod(sv[17],0);
 	fi_rbc = stod(sv[18],0);
 	fi_pla = stod(sv[19],0);

 	//double r = tao_co2_rbc_pla/tao_co2_pla_rbc;
 	//cout<<r<<endl;
 	//r = tao_hco3_rbc_pla/tao_hco3_pla_rbc;
 	//cout<<r<<endl;

}

//--------------------------------------------------------------
void solver_lumped::init_lum_tissueCO2(){
    tissueCO2v.clear();
    tissueCO2_save.clear();

    int t;
    //only one pulmonary capillary is allowed in a lumped model
	for(int i=0;i<HBsatlum ->D0_edges.size();i++){
		if(HBsatlum ->D0_edges[i]->is_per_capillary){
			t=HBsatlum ->D0_edges[i]->nx;
			//CO2 transport initialization
    		tissueCO2s = init_tissueCO2;
    		tissueCO2v.assign( t , init_tissueCO2);
		}
	}
}


//--------------------------------------------------------------
void D0_transport::prescribe_node_fi_CO2(TransportType TType, double& finode){
	switch(TType){

	case CO2_pla:
		finode = 0.02635920; //m3 O2/ m3 pla
		break;

	case CO2_rbc:
		finode = 0.02986472; //m3 O2/ m3 rbc cytoplasm
		break;

	case HCO3_pla:
		finode = 0.73471517; //m3 O2/ m3 pla
		break;

	case HCO3_rbc:
		finode = 0.08811732; //m3 O2/ m3 rbc cytoplasm
		break;

	case HbCO2:
		finode = 0.04748449; //m3 O2/ m3 rbc cytoplasm
		break;
	}

}


//--------------------------------------------------------------------------------------------------
void solver_lumped::capillary_CO2_transport(double dt){

	if(do_pul_O2_rtansport){

		pulmonary_CO2transport( dt);
	}

    if(do_per_O2_rtansport){

    	CO2transport( dt);
    }

};


//--------------------------------------------------------------------------------------------------
void solver_lumped::pulmonary_CO2transport(double dt){
	//only one capillary is allowed per 0D model.

	vector<double> CO2_pla_old = pul_cap_CO2_pla->fi;
	vector<double>& CO2_pla = pul_cap_CO2_pla->fi;

	vector<double> CO2_rbc_old = pul_cap_CO2_rbc->fi;
	vector<double>& CO2_rbc = pul_cap_CO2_rbc->fi;

	vector<double> HCO3_pla_old = pul_cap_HCO3_pla->fi;
	vector<double>& HCO3_pla = pul_cap_HCO3_pla->fi;

	vector<double> HCO3_rbc_old = pul_cap_HCO3_rbc->fi;
	vector<double>& HCO3_rbc = pul_cap_HCO3_rbc->fi;

	vector<double> HbCO2_old = pul_cap_HbCO2->fi;
	vector<double>& HbCO2 = pul_cap_HbCO2->fi;


	//double r;
	//cin>>r;
	//vfr_edge, A, nx are the same for these D0_edges
	int n = pul_cap_CO2_pla-> nx;
	vector<double> K_pul_v = sin_2(K_pul_scale_CO2, n);
	double dx = pul_cap_CO2_pla-> dx;
	double v = pul_cap_CO2_pla->corr_edge->vfr/pul_cap_CO2_pla->A*ml_to_m3;

	//double r;
	//cin>>r;
	//BCs
	double CO2_pla_n_s = pul_cap_CO2_pla->node_start->CO2_pla_n; //node start
	double CO2_pla_n_e = pul_cap_CO2_pla->node_end->CO2_pla_n; //node end
	double CO2_rbc_n_s = pul_cap_CO2_rbc->node_start->CO2_rbc_n;
	double CO2_rbc_n_e = pul_cap_CO2_rbc->node_end->CO2_rbc_n;
	double HCO3_pla_n_s = pul_cap_HCO3_pla->node_start->HCO3_pla_n;
	double HCO3_pla_n_e = pul_cap_HCO3_pla->node_end->HCO3_pla_n;
	double HCO3_rbc_n_s = pul_cap_HCO3_rbc->node_start->HCO3_rbc_n;
	double HCO3_rbc_n_e = pul_cap_HCO3_rbc->node_end->HCO3_rbc_n;
	double HbCO2_n_s = pul_cap_HbCO2->node_start->HbCO2_n;
	double HbCO2_n_e = pul_cap_HbCO2->node_end->HbCO2_n;


	//capillary
    for (int i = 1; i < n - 1; i++) {
        double CO2_pla_der;
        double CO2_rbc_der;
        double HCO3_pla_der;
        double HCO3_rbc_der;
        double HbCO2_der;

        if (v > 0.) {
            CO2_pla_der = (CO2_pla_old[i] - CO2_pla_old[i-1])/dx;
            CO2_rbc_der = (CO2_rbc_old[i] - CO2_rbc_old[i-1])/dx;
            HCO3_pla_der = (HCO3_pla_old[i] - HCO3_pla_old[i-1])/dx;
            HCO3_rbc_der = (HCO3_rbc_old[i] - HCO3_rbc_old[i-1])/dx;
            HbCO2_der = (HbCO2_old[i] - HbCO2_old[i-1])/dx;
        }
        else {
            CO2_pla_der = (CO2_pla_old[i+1] - CO2_pla_old[i])/dx;
            CO2_rbc_der = (CO2_rbc_old[i+1] - CO2_rbc_old[i])/dx;
            HCO3_pla_der = (HCO3_pla_old[i+1] - HCO3_pla_old[i])/dx;
            HCO3_rbc_der = (HCO3_rbc_old[i+1] - HCO3_rbc_old[i])/dx;
            HbCO2_der = (HbCO2_old[i+1] - HbCO2_old[i])/dx;
        }


        double Eta_hco3 = eta_hco3(HCO3_rbc_old[i], HCO3_pla_old[i], CO2_pla_old[i], CO2_rbc_old[i]);
        double Eta_hb = eta_hb( HbCO2_old[i], CO2_pla_old[i], CO2_rbc_old[i]);



        //CO2_pla[i] = CO2_pla_old[i] - dt/(fi_pla*fi_c)*v*CO2_pla_der*ksi_c*ksi_pla - dt/tao_co2_pla_rbc*(CO2_pla_old[i] - CO2_rbc_old[i]*alpha_co2_pla/alpha_co2_rbc) + dt*K_pul_scale_CO2*(PCO2_alveolar - CO2_pla_old[i]/alpha_co2_pla);
        CO2_pla[i] = CO2_pla_old[i] - dt/(fi_pla*fi_c)*v*CO2_pla_der*ksi_c*ksi_pla - dt/tao_co2_pla_rbc*(CO2_pla_old[i] - CO2_rbc_old[i]*alpha_co2_pla/alpha_co2_rbc) + dt*K_pul_v[i]*(PCO2_alveolar - CO2_pla_old[i]/alpha_co2_pla);

        CO2_rbc[i] = CO2_rbc_old[i] - dt/fi_rbc*v*CO2_rbc_der*ksi_rbc + dt/tao_co2_rbc_pla*(CO2_pla_old[i]*alpha_co2_rbc/alpha_co2_pla - CO2_rbc_old[i]) - dt/tao_hco3*Eta_hco3 - dt/tao_hbco2*Eta_hb;
        HCO3_pla[i] = HCO3_pla_old[i] - dt/fi_pla*v*HCO3_pla_der*ksi_pla - dt/tao_hco3_pla_rbc*(HCO3_pla_old[i]-HCO3_rbc_old[i]/alpha_hco3_rbc*alpha_hco3_pla);
        HCO3_rbc[i] = HCO3_rbc_old[i] - dt/fi_rbc*v*HCO3_rbc_der*ksi_rbc + dt/tao_hco3_rbc_pla*(HCO3_pla_old[i]*alpha_hco3_rbc/alpha_hco3_pla-HCO3_rbc_old[i]) + dt/tao_hco3*Eta_hco3;
        HbCO2[i] = HbCO2_old[i] - v*dt*HbCO2_der + dt/tao_hbco2*Eta_hb;
    }

    //BCs
    if (v > 0.) {
        double CO2_pla_der = (CO2_pla_old[n - 1] - CO2_pla_old[n - 2])/dx;
        double CO2_rbc_der = (CO2_rbc_old[n - 1] - CO2_rbc_old[n - 2])/dx;
        double HCO3_pla_der = (HCO3_pla_old[n - 1] - HCO3_pla_old[n - 2])/dx;
        double HCO3_rbc_der = (HCO3_rbc_old[n - 1] - HCO3_rbc_old[n - 2])/dx;
        double HbCO2_der = (HbCO2_old[n - 1] - HbCO2_old[n - 2])/dx;

        double Eta_hco3 = eta_hco3(HCO3_rbc_old[n - 1], HCO3_pla_old[n - 1], CO2_pla_old[n - 1], CO2_rbc_old[n - 1]);
        double Eta_hb = eta_hb( HbCO2_old[n - 1], CO2_pla_old[n - 1] , CO2_rbc_old[n - 1]);


        //CO2_pla[n-1] = CO2_pla_old[n-1] - dt/(fi_pla*fi_c)*v*CO2_pla_der*ksi_c*ksi_pla - dt/tao_co2_pla_rbc*(CO2_pla_old[n-1] - CO2_rbc_old[n-1]*alpha_co2_pla/alpha_co2_rbc) + dt*K_pul_scale_CO2*(PCO2_alveolar - CO2_pla_old[n-1]/alpha_co2_pla);
        CO2_pla[n-1] = CO2_pla_old[n-1] - dt/(fi_pla*fi_c)*v*CO2_pla_der*ksi_c*ksi_pla - dt/tao_co2_pla_rbc*(CO2_pla_old[n-1] - CO2_rbc_old[n-1]*alpha_co2_pla/alpha_co2_rbc) + dt*K_pul_v[n-1]*(PCO2_alveolar - CO2_pla_old[n-1]/alpha_co2_pla);
        
        CO2_rbc[n-1] = CO2_rbc_old[n-1] - dt/fi_rbc*v*CO2_rbc_der*ksi_rbc + dt/tao_co2_rbc_pla*(CO2_pla_old[n-1]*alpha_co2_rbc/alpha_co2_pla - CO2_rbc_old[n-1]) - dt/tao_hco3*Eta_hco3 - dt/tao_hbco2*Eta_hb;
        HCO3_pla[n-1] = HCO3_pla_old[n-1] - dt/fi_pla*v*HCO3_pla_der*ksi_pla - dt/tao_hco3_pla_rbc*(HCO3_pla_old[n-1]-HCO3_rbc_old[n-1]/alpha_hco3_rbc*alpha_hco3_pla);
        HCO3_rbc[n-1] = HCO3_rbc_old[n-1] - dt/fi_rbc*v*HCO3_rbc_der*ksi_rbc + dt/tao_hco3_rbc_pla*(HCO3_pla_old[n-1]*alpha_hco3_rbc/alpha_hco3_pla-HCO3_rbc_old[n-1]) + dt/tao_hco3*Eta_hco3;
        HbCO2[n-1] = HbCO2_old[n-1] - v*dt*HbCO2_der+dt/tao_hbco2*Eta_hb;

        CO2_pla[0] = CO2_pla_n_s;
        CO2_rbc[0] = CO2_rbc_n_s;
        HCO3_pla[0] = HCO3_pla_n_s;
        HCO3_rbc[0] = HCO3_rbc_n_s;
        HbCO2[0] = HbCO2_n_s;


    }
    else {
        double CO2_pla_der = (CO2_pla_old[1] - CO2_pla_old[0])/dx;
        double CO2_rbc_der = (CO2_rbc_old[1] - CO2_rbc_old[0])/dx;
        double HCO3_pla_der = (HCO3_pla_old[1] - HCO3_pla_old[0])/dx;
        double HCO3_rbc_der = (HCO3_rbc_old[1] - HCO3_rbc_old[0])/dx;
        double HbCO2_der = (HbCO2_old[1] - HbCO2_old[0])/dx;

        double Eta_hco3 = eta_hco3(HCO3_rbc_old[0], HCO3_pla_old[0], CO2_pla_old[0], CO2_rbc_old[0]);
        double Eta_hb = eta_hb( HbCO2_old[0], CO2_pla_old[0], CO2_rbc_old[0]);
        

        //CO2_pla[0] = CO2_pla_old[0] - dt/(fi_pla*fi_c)*v*CO2_pla_der*ksi_c*ksi_pla - dt/tao_co2_pla_rbc*(CO2_pla_old[0] - CO2_rbc_old[0]*alpha_co2_pla/alpha_co2_rbc) + dt*K_pul_scale_CO2*(PCO2_alveolar - CO2_pla_old[0]/alpha_co2_pla);
        CO2_pla[0] = CO2_pla_old[0] - dt/(fi_pla*fi_c)*v*CO2_pla_der*ksi_c*ksi_pla - dt/tao_co2_pla_rbc*(CO2_pla_old[0] - CO2_rbc_old[0]*alpha_co2_pla/alpha_co2_rbc) + dt*K_pul_v[0]*(PCO2_alveolar - CO2_pla_old[0]/alpha_co2_pla);
        
        CO2_rbc[0] = CO2_rbc_old[0] - dt/fi_rbc*v*CO2_rbc_der*ksi_rbc + dt/tao_co2_rbc_pla*(CO2_pla_old[0]*alpha_co2_rbc/alpha_co2_pla - CO2_rbc_old[0]) - dt/tao_hco3*Eta_hco3 - dt/tao_hbco2*Eta_hb;
        HCO3_pla[0] = HCO3_pla_old[0] - dt/fi_pla*v*HCO3_pla_der*ksi_pla - dt/tao_hco3_pla_rbc*(HCO3_pla_old[0]-HCO3_rbc_old[0]/alpha_hco3_rbc*alpha_hco3_pla);
        HCO3_rbc[0] = HCO3_rbc_old[0] - dt/fi_rbc*v*HCO3_rbc_der*ksi_rbc + dt/tao_hco3_rbc_pla*(HCO3_pla_old[0]*alpha_hco3_rbc/alpha_hco3_pla-HCO3_rbc_old[0]) + dt/tao_hco3*Eta_hco3;
        HbCO2[0] = HbCO2_old[0] - v*dt*HbCO2_der+dt/tao_hbco2*Eta_hb;

        CO2_pla[n-1] = CO2_pla_n_e;
        CO2_rbc[n-1] = CO2_rbc_n_e;
        HCO3_pla[n-1] = HCO3_pla_n_e;
        HCO3_rbc[n-1] = HCO3_rbc_n_e;
        HbCO2[n-1] = HbCO2_n_e;

    }




};

//-------------------------------------------------------------------------------------
struct TransportConfig {
    std::string edge_in_list_name;
    std::string edge_out_list_name;
    std::string node_fi_property_name;
};

//-------------------------------------------------------------------------------------
std::unordered_map<TransportType, TransportConfig> transport_configs = {
    {RBC, {"D0_edges_in_RBC", "D0_edges_out_RBC", "RBC_fi0Dn"}},
    {HB_O2_saturation, {"D0_edges_in_HBsat", "D0_edges_out_HBsat", "HBsat_0Dn"}},
    {C_Plasma_O2, {"D0_edges_in_PlasmaO2", "D0_edges_out_PlasmaO2", "PlasmaO2_0Dn"}},
    {CO2_pla, {"D0_edges_in_CO2_pla", "D0_edges_out_CO2_pla", "CO2_pla_n"}},
    {CO2_rbc, {"D0_edges_in_CO2_rbc", "D0_edges_out_CO2_rbc", "CO2_rbc_n"}},
    {HCO3_pla, {"D0_edges_in_HCO3_pla", "D0_edges_out_HCO3_pla", "HCO3_pla_n"}},
    {HCO3_rbc, {"D0_edges_in_HCO3_rbc", "D0_edges_out_HCO3_rbc", "HCO3_rbc_n"}},
    {HbCO2, {"D0_edges_in_HbCO2", "D0_edges_out_HbCO2",  "HbCO2_n"}}
};

//-------------------------------------------------------------------------------------
// Helper function to get the edge list from a node
std::vector<D0_edge*>& get_edge_list(solver_lumped::node* node, const std::string& list_name) {
    if (list_name == "D0_edges_in_RBC") return node->D0_edges_in_RBC;
    if (list_name == "D0_edges_out_RBC") return node->D0_edges_out_RBC;
    if (list_name == "D0_edges_in_HBsat") return node->D0_edges_in_HBsat;
    if (list_name == "D0_edges_out_HBsat") return node->D0_edges_out_HBsat;
    if (list_name == "D0_edges_in_PlasmaO2") return node->D0_edges_in_PlasmaO2;
    if (list_name == "D0_edges_out_PlasmaO2") return node->D0_edges_out_PlasmaO2;
    if (list_name == "D0_edges_in_CO2_pla") return node->D0_edges_in_CO2_pla;
    if (list_name == "D0_edges_out_CO2_pla") return node->D0_edges_out_CO2_pla;
    if (list_name == "D0_edges_in_CO2_rbc") return node->D0_edges_in_CO2_rbc;
    if (list_name == "D0_edges_out_CO2_rbc") return node->D0_edges_out_CO2_rbc;
    if (list_name == "D0_edges_in_HCO3_pla") return node->D0_edges_in_HCO3_pla;
    if (list_name == "D0_edges_out_HCO3_pla") return node->D0_edges_out_HCO3_pla;
    if (list_name == "D0_edges_in_HCO3_rbc") return node->D0_edges_in_HCO3_rbc;
    if (list_name == "D0_edges_out_HCO3_rbc") return node->D0_edges_out_HCO3_rbc;
    if (list_name == "D0_edges_in_HbCO2") return node->D0_edges_in_HbCO2;
    if (list_name == "D0_edges_out_HbCO2") return node->D0_edges_out_HbCO2;
    throw std::runtime_error("Unknown edge list: " + list_name);
}

//-------------------------------------------------------------------------------------
// Helper function to get the node property
double& get_node_fi(solver_lumped::node* node, const std::string& property_name) {
    if (property_name == "RBC_fi0Dn") return node->RBC_fi0Dn;
    if (property_name == "HBsat_0Dn") return node->HBsat_0Dn;
    if (property_name == "PlasmaO2_0Dn") return node->PlasmaO2_0Dn;
    if (property_name == "CO2_pla_n") return node->CO2_pla_n;
    if (property_name == "CO2_rbc_n") return node->CO2_rbc_n;
    if (property_name == "HCO3_pla_n") return node->HCO3_pla_n;
    if (property_name == "HCO3_rbc_n") return node->HCO3_rbc_n;
    if (property_name == "HbCO2_n") return node->HbCO2_n;

    throw std::runtime_error("Unknown node property: " + property_name);
}

//-------------------------------------------------------------------------------------
void D0_transport::update_nodes(solver_lumped& lum_mod) {
    // Get the configuration for the current transport type
    auto config_it = transport_configs.find(TType);
    if (config_it == transport_configs.end()) {
        std::cerr << "Error: Unknown transport type!" << std::endl;
        return;
    }

    const auto& config = config_it->second;

    for (int i = 0; i < lum_mod.nodes.size(); i++) {
        if (lum_mod.nodes[i]->is_master_node) continue; // Skip master nodes

        double q = 0.0;
        double c = 0.0; // concentration

        // Incoming edges
        auto& edges_in = get_edge_list(lum_mod.nodes[i], config.edge_in_list_name);
        for (auto* edge : edges_in) {
            double Q = edge->corr_edge->vfr * edge->corr_edge->is_open;
            if (Q > 0.0) {
                q += Q;
                c += Q * edge->fi.back();
            }
        }

        // Outgoing edges
        auto& edges_out = get_edge_list(lum_mod.nodes[i], config.edge_out_list_name);
        for (auto* edge : edges_out) {
            double Q = edge->corr_edge->vfr * edge->corr_edge->is_open;
            if (Q < 0.0) {
                q -= Q;
                c -= Q * edge->fi[0];
            }
        }

        // Update node concentration
        if (q != 0.0) {
            get_node_fi(lum_mod.nodes[i], config.node_fi_property_name) = c / q;
        }
    }
}

//-------------------------------------------------------------------------------------
void D0_transport::connect_0D_edges(solver_lumped& lum_mod) {
    for (int i = 0; i < D0_edges.size(); i++) {
        std::string ns = D0_edges[i]->node_s_name;
        std::string ne = D0_edges[i]->node_e_name;

        int index_start = -1;
        int index_end = -1;

        for (int j = 0; j < lum_mod.nodes.size(); j++) {
            if (ns == lum_mod.nodes[j]->name) {
                index_start = j;
                // Get the configuration for the current transport type
                auto config_it = transport_configs.find(TType);
                if (config_it == transport_configs.end()) {
                    std::cerr << "Error: Unknown transport type!" << std::endl;
                    exit(-1);
                }
                const auto& config = config_it->second;
                // Add the edge to the outgoing edge list of the start node
                get_edge_list(lum_mod.nodes[j], config.edge_out_list_name).push_back(D0_edges[i]);
            }

            if (ne == lum_mod.nodes[j]->name) {
                index_end = j;
                // Get the configuration for the current transport type
                auto config_it = transport_configs.find(TType);
                if (config_it == transport_configs.end()) {
                    std::cerr << "Error: Unknown transport type!" << std::endl;
                    exit(-1);
                }
                const auto& config = config_it->second;
                // Add the edge to the incoming edge list of the end node
                get_edge_list(lum_mod.nodes[j], config.edge_in_list_name).push_back(D0_edges[i]);
            }
        }

        if (index_start < 0 || index_end < 0) {
            std::cout << "Transport node " << ns << " or " << ne << " does not exist in the lumped model." << std::endl;
            exit(-1);
        }

        D0_edges[i]->node_start = lum_mod.nodes[index_start];
        D0_edges[i]->node_end = lum_mod.nodes[index_end];
    }
}

const std::string& D0_transport::get_node_fi_property_name(TransportType type) {
    auto config_it = transport_configs.find(type);
    if (config_it == transport_configs.end()) {
        throw std::runtime_error("Unknown transport type!");
    }
    return config_it->second.node_fi_property_name;
}

//-------------------------------------------------------------------------------------
void D0_edge::update_diode() {
    if (corr_edge->is_open) {
        try {
            const std::string& property_name = D0_transport::get_node_fi_property_name(TType);
            fi[0] = get_node_fi(node_start, property_name);
        } catch (const std::runtime_error& e) {
            std::cerr << "Error in update_diode: " << e.what() << std::endl;
        }
    }
}


//-------------------------------------------------------------------------------------
void D0_edge::virt1D(double dt) {
    double v = corr_edge->vfr / A * ml_to_m3;
    std::vector<double> fi_old = fi;

    // Update interior points
    for (int i = 1; i < nx - 1; i++) {
        if (v > 0.) {
            fi[i] = fi_old[i] - v * dt / dx * (fi_old[i] - fi_old[i - 1]);
        }
        else {
            fi[i] = fi_old[i] - v * dt / dx * (fi_old[i + 1] - fi_old[i]);
        }
    }

    // Boundary conditions
    if (v > 0.) {
        fi[nx - 1] = fi_old[nx - 1] - v * dt / dx * (fi_old[nx - 1] - fi_old[nx - 2]);

        try {
            const std::string& property_name = D0_transport::get_node_fi_property_name(TType);
            fi[0] = get_node_fi(node_start, property_name);
        } catch (const std::runtime_error& e) {
            std::cerr << "Error in virt1D: " << e.what() << std::endl;
        }

    }
    else {
        try {
            const std::string& property_name = D0_transport::get_node_fi_property_name(TType);
            fi[nx - 1] = get_node_fi(node_end, property_name);
        } catch (const std::runtime_error& e) {
            std::cerr << "Error in virt1D: " << e.what() << std::endl;
        }

        fi[0] = fi_old[0] - v * dt / dx * (fi_old[1] - fi_old[0]);

    }
}



//Co2 control in brain 

//--------------------------------------------------------------
void solver_lumped::CO2_response(double t_act)
{

	// time step
	double dt = t_act - time.back();

	double P_CO2 = P_CO2_ave->average.back(); //partial pressure of arterial co2 locally

	// actuator signal
	//cout<<P_CO2<<"  "<<CO2_ref<<endl<<endl;
	cout<<t_act<<endl;
	x_CO2 = x_CO2 + dt / tao_CO2 * (- x_CO2 + G_CO2 * (P_CO2 - CO2_ref)/CO2_ref); // le kell normálni

}


//--------------------------------------------------------------
void solver_lumped::load_perif_baroreflex(string filename){
	do_perif_baroreflex = true;

	ifstream file_in;
	string file_name = input_folder_path + '/' + filename + ".txt";
	file_in.open(file_name);
	string line;
	if(file_in.is_open()){
		while(getline(file_in,line)){
			line.erase(remove(line.begin(), line.end(), ' '), line.end());
			line.erase(remove(line.begin(), line.end(), '\n'), line.end());
			line.erase(remove(line.begin(), line.end(), '\r'), line.end());

			// separating strings to vector by comma
			vector<string> sv = separate_line(line);
			x_bar_time.push_back(stod(sv[0],0));
			x_bar_vect.push_back(stod(sv[1],0));

		}
	}
	else{
		cout << "! ERROR !" << endl << " File is not open when calling load_main_csv() function!!! file: " << file_name << "\n" << endl;
	}
}


void solver_lumped::perif_baroreflex()
{
    double t = time.back();

    if (t <= x_bar_time.front()) {
        x_bar = x_bar_vect.front();
        return;
    }

    if (t >= x_bar_time.back()) {
        x_bar = x_bar_vect.back();
        return;
    }

    //interval
    while (index_of_time + 1 < x_bar_time.size() &&
           t >= x_bar_time[index_of_time + 1]) {
        index_of_time++;
    }

    double t0 = x_bar_time[index_of_time];
    double t1 = x_bar_time[index_of_time + 1];
    double y0 = x_bar_vect[index_of_time];
    double y1 = x_bar_vect[index_of_time + 1];

    x_bar = y0 + (t - t0) * (y1 - y0) / (t1 - t0);
}