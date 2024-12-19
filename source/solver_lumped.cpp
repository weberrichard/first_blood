#include "solver_lumped.h"

//--------------------------------------------------------------
solver_lumped::solver_lumped(string a_name, string a_folder)
{
	name = a_name;
	input_folder_path = a_folder;
}

solver_lumped::~solver_lumped(){}

//--------------------------------------------------------------
void solver_lumped::initialization(double hr)
{
	// setting sizes
	number_of_nodes = nodes.size();
	number_of_edges = edges.size();
	number_of_master = boundary_model_node.size();

	// clearing master boundary indices
	boundary_indices.clear();

	// heart rate
	heart_rate = hr; // from Charlton2019
	time_period = 60./heart_rate;

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
	q_ave = new time_average();
	p_ave = new time_average();
	C_ave = new time_average();
	R_fact = new time_average();
	x_myo_ts = new time_average();

	//for metabolic response
	Ct_ave = new time_average();
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

		// updating parameters: applying control effects
	if(t_act>3.*time_period)
	{
		update_parameters(t_act);
	}


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
		double vn = edges[0]->vfr;
		q_ave->update(tn,vn,time_period);

		vn = nodes[5]->p;
		p_ave->update(tn,vn,time_period);
		
		vn = edges[5]->par_non_SI[0];
		C_ave->update(tn,vn,time_period);

		vn = edges[0]->parameter_factor;
		R_fact->update(tn,vn,time_period);
		x_myo_ts->update(tn,x_myo,time_period);
	}


	//saving tissue O2 concentration
	if(do_lum_PlasmaO2_transport&&do_lum_HB_sat_transport&&do_lum_RBC_transport){
        tissueO2_save.push_back(tissueO2s);
    }

    //update for metabolic response
    if(do_metabolic_res){
    	double tn = time.back();
    	Ct_ave->update(tn, tissueO2s, time_period);
    }

}


//--------------------------------------------------------------
void solver_lumped::update_parameters(double t_act)
{
	autoregulation(t_act);
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
		cout << "\n!!!WARNING!!!\n solver_lumped::edge_id_to_index function\n Node is not existing, edge_id: " << edge_id << "\n Continouing..." << endl;
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
{	
	// normalized version
	double tn = t * heart_rate/60.;

	// making the elastance periodic
	while(tn>1.)
	{
		tn -= 1.;
	}

	double En = 17.4073 * pow(tn,1.9) / (1.+11.2305*pow(tn,1.9)) * 1. / (1.+1.6658e7*pow(tn,21.9));

	//double En = 1.55*pow(tn/0.7,1.9)/(1.+pow(tn/0.7,1.9)) * (1./(1.+pow(tn/1.17,21.9)));

	double E = (par[0]-par[1])*En + par[1];

	// E = E*mmHg_to_Pa*1.e6; // mmHg/ml to SI: Pa/m3 

	return E;
}

//--------------------------------------------------------------
double solver_lumped::elastance_derived(double t, vector<double> par)
{
	// normalized version
	double tn = t * heart_rate/60.;

	// making the elastance periodic
	while(tn>1.)
	{
		tn -= 1.;
	}

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
double solver_lumped::delta_V(int edge_index, int node_index){
	double C_ref = edges[edge_index]->parameter[0]; //everything is in SI
	double dp = nodes[node_index]->p*mmHg_to_Pa - atmospheric_pressure;
	return C_ref * dp;
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
void D0_transport::update_fi(double dt,solver_lumped& lum_mod){
	update_nodes( lum_mod);
	update_edges( dt);
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

//	if(name=="p10"){
//	Mmax = 2.7e-4*2.0;
//}

	//vfr_edge, A, nx are the same for these D0_edges
	int n = per_cap_RBC-> nx;
	double dx = per_cap_RBC-> dx;
	double v = per_cap_RBC->vfr_edge->vfr/per_cap_RBC->A*ml_to_m3;

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
	double v = pul_cap_BH->vfr_edge->vfr/pul_cap_BH->A * ml_to_m3;
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

    //only one pulmonary capillary is allowed in a lumped model but noth both
    int t;
	for(int i=0;i<HBsatlum ->D0_edges.size();i++){
		if(HBsatlum ->D0_edges[i]->is_pul_capillary){
			t=HBsatlum ->D0_edges[i]->nx;
		}
	}

	for(int i=0;i<HBsatlum ->D0_edges.size();i++){
		if(HBsatlum ->D0_edges[i]->is_per_capillary){
			t=HBsatlum ->D0_edges[i]->nx;
		}
	}

    //O2 transport initialization
    tissueO2s = init_tissueO2;
    tissueO2v.assign( t , init_tissueO2);
}


//--------------------------------------------------------------
double solver_lumped::turn_source(double t){
	return 1/(1 + exp(-(t-10.)));
}


//--------------------------------------------------------------
void solver_lumped::save_tissueO2(string folder_name, const vector<double>& st, const vector<double>& time){
	   if (do_lum_PlasmaO2_transport&&do_lum_HB_sat_transport&&do_lum_RBC_transport){

		string file_name = folder_name;
		HBsatlum->save_vector(file_name, tissueO2_save, time);
   }
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
    Mmax = stod(sv[8],0);
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
	K_pul_O2 = stod(sv[2],0);
	taoO2_p = stod(sv[3],0);
	K_pul_scale =stod(sv[4],0) ;
}


//--------------------------------------------------------------
void solver_lumped::metabolic_response(double t_act)
{

	// time step
	double dt = t_act - time.back();

	double Ct = Ct_ave->average.back();//p_ave->average.back();

	// actuator signal
	x_met = x_met + dt / tao_met * (- x_met + G_met * (Ct - Ct_ref)/Ct_ref);

}


double solver_lumped::vessel_dilation(int edgeindex){
	//return 1.;
	return pow( edges[edgeindex]->parameter_factor ,-0.5);
}


//--------------------------------------------------------------------------------------------------
D0_edge::D0_edge(string D0_name, double L, double A, int  nx, TransportType TType, double init, string node_s_name, string node_e_name, string diode_name, string vfr_edge_name):D0_name(D0_name),L(L),A(A),nx(nx),TType(TType),init(init),node_s_name(node_s_name),node_e_name(node_e_name),vfr_edge_name(vfr_edge_name),diode_name(diode_name){

	dx = L/(nx-1);
	fi_start.clear();
	fi_end.clear();
	fi.clear();

	fi.assign(nx, init);
	if(do_save_memory){
		save();}

}


//--------------------------------------------------------------------------------------------------
void D0_transport::update_nodes(solver_lumped& lum_mod){
	for(int i=0; i<lum_mod.nodes.size() ;i++){
		double q=0.;
		double c=0.; //concantration

		if(!lum_mod.nodes[i]->is_master_node){ //master nodes are handled separately in a different function
			switch(TType){
				case RBC:
					//incoming edges
					for(int j=0; j< lum_mod.nodes[i]->D0_edges_in_RBC.size() ; j++ ){
						D0_edge* d = lum_mod.nodes[i]->D0_edges_in_RBC[j];
						double Q = d->vfr_edge->vfr;
						if(Q > 0.){
							q += Q;
							c += Q * d->fi.back();
						}
					}

					//outgoing edges
					for(int j=0; j< lum_mod.nodes[i]->D0_edges_out_RBC.size() ; j++ ){
						D0_edge* d = lum_mod.nodes[i]->D0_edges_out_RBC[j];
						double Q = d->vfr_edge->vfr;
						if(Q < 0.){
							q -= Q;
							c -= Q * d->fi[0];
						}
					}
					if(q !=0. ){lum_mod.nodes[i]->RBC_fi0Dn = c/q;}
				break;

				case C_Plasma_O2:
					//incoming edges
					for(int j=0; j< lum_mod.nodes[i]->D0_edges_in_PlasmaO2.size() ; j++ ){
						D0_edge* d = lum_mod.nodes[i]->D0_edges_in_PlasmaO2[j];
						double Q = d->vfr_edge->vfr;
						if(Q > 0.){
							q += Q;
							c += Q * d->fi.back();
						}
					}

					//outgoing edges
					for(int j=0; j< lum_mod.nodes[i]->D0_edges_out_PlasmaO2.size() ; j++ ){
						D0_edge* d = lum_mod.nodes[i]->D0_edges_out_PlasmaO2[j];
						double Q = d->vfr_edge->vfr;
						if(Q < 0.){
							q -= Q;
							c -= Q * d->fi[0];
						}
					}
					if(q !=0. ){lum_mod.nodes[i]->PlasmaO2_0Dn = c/q;}
				break;

				case HB_O2_saturation:
					//incoming edges
					for(int j=0; j< lum_mod.nodes[i]->D0_edges_in_HBsat.size() ; j++ ){
						D0_edge* d = lum_mod.nodes[i]->D0_edges_in_HBsat[j];
						double Q = d->vfr_edge->vfr;
						if(Q > 0.){
							q += Q;
							c += Q * d->fi.back();
						}
					}

					//outgoing edges
					for(int j=0; j< lum_mod.nodes[i]->D0_edges_out_HBsat.size() ; j++ ){
						D0_edge* d = lum_mod.nodes[i]->D0_edges_out_HBsat[j];
						double Q = d->vfr_edge->vfr;
						if(Q < 0.){
							q -= Q;
							c -= Q * d->fi[0];
						}
					}
					if(q !=0. ){lum_mod.nodes[i]->HBsat_0Dn = c/q;}
				break;
			}
		}
		//capacitances modelling dilation and contraction
	}
}


//--------------------------------------------------------------------------------------------------
void D0_transport::update_edges( double dt){
//capillary edges are updates from solver_lumped
	for(int i=0; i< D0_edges.size(); i++ ){

		if((!D0_edges[i]->is_per_capillary && !D0_edges[i]->is_pul_capillary) || !do_tissue_transport){
			if(D0_edges[i]->is_diode){
				D0_edges[i]->update_diode();
			}
			else{
				D0_edges[i]->virt1D(dt);
			}
		}
	}
}


//--------------------------------------------------------------------------------------------------
void D0_edge::virt1D(double dt){
	double v = vfr_edge->vfr/A * ml_to_m3;
	
	vector<double> fi_old = fi;

    for (int i = 1; i < nx - 1; i++) {

        if (v > 0.) {
            fi[i] = fi_old[i] - v * dt / dx * (fi_old[i] - fi_old[i - 1]);
        }
        else {
            fi[i] = fi_old[i] - v * dt / dx * (fi_old[i + 1] - fi_old[i]);
        }
    }

    //BC
    
    if (v > 0.) {

        fi[nx - 1] = fi_old[nx - 1] - v * dt / dx * (fi_old[nx - 1] - fi_old[nx - 2]);

        double fiStartNode;
        switch(TType){
    case RBC:
    	fiStartNode = node_start->RBC_fi0Dn;
    	break;

    case C_Plasma_O2:
    	fiStartNode = node_start->PlasmaO2_0Dn;
    	break;

    case HB_O2_saturation:
    	fiStartNode = node_start->HBsat_0Dn;
    	break;}


        fi[0] = fiStartNode;
    }
    else {

    	double fiEndNode;
        switch(TType){
    case RBC:
    	fiEndNode = node_end->RBC_fi0Dn;
    	break;

    case C_Plasma_O2:
    	fiEndNode = node_end->PlasmaO2_0Dn;
    	break;

    case HB_O2_saturation:
    	fiEndNode = node_end->HBsat_0Dn;
    	break;}

        fi[nx - 1] = fiEndNode; 
        fi[0] = fi_old[0] - v * dt / dx * (fi_old[1] - fi_old[0]);
    }

}


void D0_edge::save(){

	fi_start.push_back(fi[0]);
	fi_end.push_back(fi.back());
}

void D0_edge::update_diode(){

	if(D0_diode->is_open){
		switch(TType){
		case RBC:
		fi[1] = node_start->RBC_fi0Dn;
		fi[0] = fi[1];
		node_end->RBC_fi0Dn = fi[1];
		break;

		case HB_O2_saturation:
		fi[1] = node_start->HBsat_0Dn;
		fi[0] = fi[1];
		node_end->HBsat_0Dn = fi[1];
		break;

		case C_Plasma_O2:
		fi[1] = node_start->PlasmaO2_0Dn;
		fi[0] = fi[1];
		node_end->PlasmaO2_0Dn = fi[1];
		break;
		}
	}
}


//--------------------------------------------------------------------------------------------------
void D0_transport::connect_0D_edges(solver_lumped& lum_mod){
	for(int i=0; i<D0_edges.size(); i++){
		string ns = D0_edges[i]->node_s_name;
		string ne = D0_edges[i]->node_e_name;

		int index_start = -1;
		int index_end = -1;
		for(int j=0;j<lum_mod.nodes.size();j++){
			if (ns == lum_mod.nodes[j]->name){
				index_start=j;

				switch(TType){
				case RBC:
				lum_mod.nodes[j]->D0_edges_out_RBC.push_back(D0_edges[i]);
				break;

				case HB_O2_saturation:
				lum_mod.nodes[j]->D0_edges_out_HBsat.push_back(D0_edges[i]);
				break;

				case C_Plasma_O2:
				lum_mod.nodes[j]->D0_edges_out_PlasmaO2.push_back(D0_edges[i]);
				break;
				}

			}
			if (ne == lum_mod.nodes[j]->name){
				index_end=j;

				switch(TType){
				case RBC:
				lum_mod.nodes[j]->D0_edges_in_RBC.push_back(D0_edges[i]);
				break;

				case HB_O2_saturation:
				lum_mod.nodes[j]->D0_edges_in_HBsat.push_back(D0_edges[i]);
				break;

				case C_Plasma_O2:
				lum_mod.nodes[j]->D0_edges_in_PlasmaO2.push_back(D0_edges[i]);
				break;}
			}
		}

		if(index_start<0 || index_end<0){
			cout<<"Transport node " <<ns<< " or " <<ne<<" does not exist in the lumped model."<<endl;
			exit(-1);
		}

		D0_edges[i]->node_start = lum_mod.nodes[index_start];
		D0_edges[i]->node_end = lum_mod.nodes[index_end];

	}
}


//--------------------------------------------------------------------------------------------------
void solver_lumped::set_0D_pointers(){

	if(do_lum_RBC_transport){
		for(int i=0; i<RBClum->D0_edges.size(); i++){
			if(RBClum->D0_edges[i]->is_diode){
				for(int j=0; j<edges.size();j++){
					if(RBClum->D0_edges[i]->diode_name == edges[j]->name ){
						RBClum->D0_edges[i]->D0_diode = edges[j];
						RBClum->D0_edges[i]->vfr_edge = edges[j];
					}
				}
			}

			else{
				for(int j=0; j<edges.size();j++){
					if(RBClum->D0_edges[i]->vfr_edge_name == edges[j]->name ){RBClum->D0_edges[i]->vfr_edge = edges[j];}
				}
			}
		}
	}


	if(do_lum_PlasmaO2_transport){
		for(int i=0; i<PlasmaO2lum->D0_edges.size(); i++){
			if(PlasmaO2lum->D0_edges[i]->is_diode){
				for(int j=0; j<edges.size();j++){
					if(PlasmaO2lum->D0_edges[i]->diode_name == edges[j]->name ){
						PlasmaO2lum->D0_edges[i]->D0_diode = edges[j];
						PlasmaO2lum->D0_edges[i]->vfr_edge = edges[j];
					}
				}
		
			}
			else{
				for(int j=0; j<edges.size();j++){
					if(PlasmaO2lum->D0_edges[i]->vfr_edge_name == edges[j]->name ){PlasmaO2lum->D0_edges[i]->vfr_edge = edges[j];}
				}
			}
		}
	}

	if(do_lum_HB_sat_transport){
		for(int i=0; i<HBsatlum->D0_edges.size(); i++){
			if(HBsatlum->D0_edges[i]->is_diode){
				for(int j=0; j<edges.size();j++){
					if(HBsatlum->D0_edges[i]->diode_name == edges[j]->name ){
						HBsatlum->D0_edges[i]->D0_diode = edges[j];
						HBsatlum->D0_edges[i]->vfr_edge = edges[j];
					}
				}
		
			}
			else{
				for(int j=0; j<edges.size();j++){
					if(HBsatlum->D0_edges[i]->vfr_edge_name == edges[j]->name ){HBsatlum->D0_edges[i]->vfr_edge = edges[j];}
				}
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
	if(do_myogenic)
	{
		//updates x_myo
		myogenic_control(t_act);
	}

	if(do_metabolic_res){
		//updates x_met
		metabolic_response(t_act);
	}

	//updates the parameter factor of the resistance
	update_R_fact();
}


//--------------------------------------------------------------------------------------------------
void solver_lumped::update_R_fact(){

	vector<int> Ridx{0,1}; // which resistors are we modifying

	double FF; 
	for(int i=0; i<Ridx.size(); i++)
	{
		double Rmax, Rmin; // calculated from r_min, r_max from "Regulation of Coronary Microvascular Resistance in Health and Disease" pic 12.2
		if((x_met + x_myo) < 0){
			//the sigmoid curve is the same for the two responses
			Rmin = sat1_met;
			Rmax = 2. - Rmin;
		}
		else{
			Rmax = sat2_met;
			Rmin = 2. - sat2_met;
		}

		double ff = 80. / (Rmax - Rmin) ;
		FF = (Rmax + Rmin * exp(- (x_met + x_myo) * ff)) / (1. + exp(-(x_met + x_myo) * ff));
		
		//cout.precision(10);
		//cout << FF << endl;

		edges[Ridx[i]]->parameter_factor = FF;
	}


}