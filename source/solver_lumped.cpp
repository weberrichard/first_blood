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
	}

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

	// setting Eigen vars
	int nm = number_of_edges + number_of_nodes + number_of_master + 2*number_of_elastance;
	A = MatrixXd::Zero(nm,nm);
	b = VectorXd::Zero(nm);

	// building model
	build_system();

	// for myogenic control
	q_ave = new time_average();
	p_ave = new time_average();
	C_ave = new time_average();
	R_fact = new time_average();
	x_myo_ts = new time_average();
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


}


//--------------------------------------------------------------
void solver_lumped::update_parameters(double t_act)
{
	if(do_myogenic)
	{
		myogenic_control(t_act);
	}
}

//--------------------------------------------------------------
void solver_lumped::myogenic_control(double t_act)
{

	// time step
	double dt = t_act - time.back();

	double p = p_ave->average.back();

	// actuator signal
	x_myo = x_myo + dt / tao * (- x_myo + G * (p - p_ref));

	vector<int> Ridx{0,1}; // which resistors are we modifying

	double FF; 
	for(int i=0; i<Ridx.size(); i++)
	{
		double Rmax, Rmin; // calculated from r_min, r_max from "Regulation of Coronary Microvascular Resistance in Health and Disease" pic 12.2
		if(x_myo < 0){
			Rmax = 1.228;
			Rmin = 0.772;
		}
		else{
			Rmax = 1.773;
			Rmin = 0.227;
		}

		double R_ref = edges[Ridx[i]]->par_non_SI[0];
		double K = (p_ref-atmospheric_pressure/mmHg_to_Pa) * (p_ref-atmospheric_pressure/mmHg_to_Pa) / R_ref;
		double ff = 10*8. * (p_ref-atmospheric_pressure/mmHg_to_Pa) / ( K * (Rmax - Rmin) * R_ref );
		FF = (Rmax + Rmin * exp(-x_myo * ff)) / (1. + exp(-x_myo * ff));
		
		//cout.precision(10);
		//cout << FF << endl;

		edges[Ridx[i]]->parameter_factor = FF;
	}

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


//-----------------------------------------
void Virt1DforLum(vector<double> &fi_old, vector<double> &fi, double v, double dt, double dx, int n, double fiStartNode, double fiEndNode) {
    vector<double> fi_tmp = fi;

    for (int i = 1; i < n - 1; i++) {

        if (v > 0.) {
            fi[i] = fi_old[i] - v * dt / dx * (fi_old[i] - fi_old[i - 1]);
        }
        else {
            fi[i] = fi_old[i] - v * dt / dx * (fi_old[i + 1] - fi_old[i]);
        }
    }

    //BC
    
    if (v > 0.) {

        fi[n - 1] = fi_old[n - 1] - v * dt / dx * (fi_old[n - 1] - fi_old[n - 2]);
        fi[0] = fiStartNode;
    }
    else {
        fi[n - 1] = fiEndNode; 
        fi[0] = fi_old[0] - v * dt / dx * (fi_old[1] - fi_old[0]);
    }

    fi_old = fi_tmp;

}


//----------------------------------------
D0_transport::D0_transport(LumpedType LType, vector<string> sv, TransportType TType):LType(LType), TType(TType) {
    switch (LType) {
    case PerifCoronary0D:
    //no idea
        break;
    case Perif0D:

    	initialization();

        L_arteriole = stod(sv[1],0);
        L_capillary = stod(sv[2],0);
        L_venulare = stod(sv[3],0);
        L_vein = stod(sv[4],0);
        A_arteriole = stod(sv[5],0);
        A_capillary = stod(sv[6],0);
        A_venulare = stod(sv[7],0);
        A_vein = stod(sv[8],0);

        //arteriole
        nx_arteriole = NX( L_arteriole, stod(sv[9],0) , 10); 
        dx_arteriole = L_arteriole / (nx_arteriole - 1);
        fi_arteriole.assign(nx_arteriole, 0.);
        fi_old_arteriole.assign(nx_arteriole, 0.);

        //capillary
        nx_capillary = NX(L_capillary, stod(sv[10],0), 5);
        dx_capillary = L_capillary / (nx_capillary - 1);
        fi_capillary.assign(nx_capillary, 0.);
        fi_old_capillary.assign(nx_capillary, 0.);

        //venulare
        nx_venulare = NX(L_venulare, stod(sv[11],0), 5);
        dx_venulare = L_venulare / (nx_venulare - 1);
        fi_venulare.assign(nx_venulare, 0.);
        fi_old_venulare.assign(nx_venulare, 0.);

        //vein
        nx_vein = NX(L_vein, stod(sv[12],0), 60);
        dx_vein = L_vein / (nx_vein - 1);
        fi_vein.assign(nx_vein, 0.);
        fi_old_vein.assign(nx_vein, 0.);

        
        save_variables();

        break;
    case Heart0D:
    	//L, A, nx
    	L_pul_art = stod(sv[1],0);
    	L_pul_vein = stod(sv[2],0);
    	A_pul_art = stod(sv[3],0);
    	A_pul_vein = stod(sv[4],0);
    	nx_pul_art = stod(sv[5],0);
    	nx_pul_vein = stod(sv[6],0);
    	dx_pul_art = L_pul_art / (nx_pul_art - 1);
    	dx_pul_vein = L_pul_vein / (nx_pul_vein - 1);

    	//nodes
    	fi_RA, fi_RV, fi_LA, fi_LV = 0.;
        fi_old_RA, fi_old_RV, fi_old_LA, fi_old_LV = 0.;
        fi_lung = 0.;

        //virtual 1D
        fi_pul_art.assign(nx_pul_art, 0.);
        fi_pul_vein.assign(nx_pul_vein, 0.);
        fi_old_pul_art.assign(nx_pul_art, 0.);
        fi_old_pul_vein.assign(nx_pul_vein, 0.);

        save_variables();

        break;
    default:
        cout << "Unknown 0D type for transport.";
    }
}


//-----------------------------------
void D0_transport::update_fi(double dt, double& masterFi, solver_lumped& lum_mod, double fi_vena_cava) {
	double AA; //changing cross-section if needed

    switch (this-> LType) {
    case PerifCoronary0D:
        //no idea
        break;
    case Perif0D:

    	switch(TType){
    	case RBC:
        
        AA = A_arteriole + lum_mod.delta_V( 5, 1)/L_arteriole;
        Virt1DforLum(fi_old_arteriole, fi_arteriole, lum_mod.edges[1]->vfr * ml_to_m3 / AA, dt, dx_arteriole, nx_arteriole, masterFi , lum_mod.nodes[2]->RBC_fi0Dn);

        AA = A_capillary + lum_mod.delta_V( 6, 2)/L_capillary;
        Virt1DforLum(fi_old_capillary, fi_capillary, lum_mod.edges[2]->vfr * ml_to_m3 / AA, dt, dx_capillary, nx_capillary, lum_mod.nodes[2]->RBC_fi0Dn, lum_mod.nodes[3]->RBC_fi0Dn);

        AA = A_venulare + lum_mod.delta_V( 7, 3)/L_venulare;
        Virt1DforLum(fi_old_venulare, fi_venulare, lum_mod.edges[3]->vfr * ml_to_m3 / AA, dt, dx_venulare, nx_venulare, lum_mod.nodes[3]->RBC_fi0Dn, lum_mod.nodes[4]->RBC_fi0Dn);

        AA = A_vein + lum_mod.delta_V( 8, 4)/L_vein;
        Virt1DforLum(fi_old_vein, fi_vein, lum_mod.edges[4]->vfr * ml_to_m3 / AA, dt, dx_vein, nx_vein, lum_mod.nodes[4]->RBC_fi0Dn, fi_vena_cava);
        break;
        }
        

        //nodes
        //not sure if the virtual 1D or the nodes should be updated first
        UpdatePerifLumNode(2, fi_arteriole.back(), fi_capillary[0], lum_mod);
        UpdatePerifLumNode(3, fi_capillary.back(), fi_venulare[0], lum_mod);
        UpdatePerifLumNode(4, fi_venulare.back(), fi_vein[0], lum_mod);

        break;
    case Heart0D:// that will be fun...

        //virt 1D for lungs...
        //pul arteries
    	AA = A_pul_art + lum_mod.delta_V( 7, 8)/ (L_pul_vein + L_pul_art);
    	if(lum_mod.edges[5]->is_open){
    		Virt1DforLum(fi_old_pul_art, fi_pul_art, lum_mod.edges[5]->vfr * ml_to_m3 / AA, dt, dx_pul_art, nx_pul_art, fi_RV, fi_lung);
        }
        else{
        	Virt1DforLum(fi_old_pul_art, fi_pul_art, lum_mod.edges[5]->vfr * ml_to_m3 / AA, dt, dx_pul_art, nx_pul_art, fi_old_pul_art[0] , fi_lung);
        }

        //pul vein
        AA = A_pul_vein + lum_mod.delta_V( 7, 8)/ (L_pul_vein + L_pul_art);
        Virt1DforLum(fi_old_pul_vein, fi_pul_vein, lum_mod.edges[6]->vfr * ml_to_m3 / AA, dt, dx_pul_vein, nx_pul_vein, fi_lung, fi_LA);

        //nodes
        //right atrium
        fi_RA = fi_vena_cava;
        //right ventricle
        if (lum_mod.edges[2]->is_open){
            fi_RV = fi_RA;
        }

        if(lum_mod.edges[5]->is_open){
        	fi_pul_art[0] = fi_RV;
        }

        //update_lung_fi(fi_lung, fi_pul_art.back(), fi_pul_vein[0], lum_mod);
        prescribe_lung_fi(lum_mod);

        //left atrium
        fi_LA = fi_pul_vein.back(); //can't flow in the other direction, must be the same
        //left ventricle
        if(lum_mod.edges[10]->is_open){
            fi_LV = fi_LA;
        }
        //master_node is upsated in TransportNodeCl::update_master_fi

        break;
    default:
        cout << "Unknown 0D type for transport.";
    }

    if(do_save_results){
        save_variables(); //not sure about it...
        }

    return;
}

//------------------------------------
//updates fi parameters of perif nodes
void D0_transport::UpdatePerifLumNode(int LumNodeIndex, double fiLeft, double fiRight, solver_lumped& lum_mod) {
    int a=0, b=0, c=0;

    //1 if q flows towards the node
    if (lum_mod.edges[LumNodeIndex]->vfr < 0.) { a = 1; } //Resistance edge right
    if (lum_mod.edges[LumNodeIndex + 7]->vfr > 0.) { b = 1; } //Inductance edge left
    if (lum_mod.edges[LumNodeIndex + 4]->vfr < 0.) { c = 1; } //Capacitance edge

    double qLeft, qRight, qDown;
    
    switch(TType){
    case RBC:

    switch (a*4 + b*2 + c) { // in case 0 nothing changes
    case 1:
        lum_mod.nodes[LumNodeIndex]->RBC_fi0Dn = fiRight;
        break;

    case 2:
        lum_mod.nodes[LumNodeIndex]->RBC_fi0Dn = fiLeft;
        break;

    case 3:
        qLeft = lum_mod.edges[LumNodeIndex + 7]->vfr;
        qDown = -1.* lum_mod.edges[LumNodeIndex + 4]->vfr;
        lum_mod.nodes[LumNodeIndex]->RBC_fi0Dn = (qLeft*fiLeft + qDown *fiRight)/(qDown + qLeft);
        break;

    case 4:
        lum_mod.nodes[LumNodeIndex]->RBC_fi0Dn = fiRight; //same az case 1
        break;

    case 5:
        lum_mod.nodes[LumNodeIndex]->RBC_fi0Dn = fiRight; //both has the same fi
        break;

    case 6:
        qLeft = lum_mod.edges[LumNodeIndex + 7]->vfr;
        qRight = -1.* lum_mod.edges[LumNodeIndex]->vfr;
        lum_mod.nodes[LumNodeIndex]->RBC_fi0Dn = (qLeft * fiLeft + qRight * fiRight) / (qRight + qLeft);
        break;
    }

    break;
    }





}

//-----------------------------------------------------------
void D0_transport::update_lung_fi(double& fi_lung, double fiLeft, double fiRight, solver_lumped& lum_mod){
	int a=0, b=0, c=0;

    //1 if q flows towards the node
    if (lum_mod.edges[5]->vfr > 0.) { a = 1; } //diode R
    if (lum_mod.edges[6]->vfr < 0.) { b = 1; } //Rp
    if (lum_mod.edges[7]->vfr > 0.) { c = 1; } //Cp

    double qLeft, qRight, qDown, fiDown;
    //cout<<a*4 + b*2 + c<<endl;

    switch (a*4 + b*2 + c) {
    case 4:
        fi_lung = fiRight;
        break;

    case 5:
        qLeft = lum_mod.edges[5]->vfr;
        qRight = lum_mod.edges[6]->vfr;//no negative sign, think it through...
        qDown = lum_mod.edges[7]->vfr;
        fiDown = (qLeft * fiLeft + qRight * fiRight) / (qRight + qLeft);

        fi_lung = (qLeft * fiLeft + qDown * fiDown) / (qDown + qLeft);
        break;

    case 6:
        qLeft = lum_mod.edges[5]->vfr;
        qRight = -1.* lum_mod.edges[6]->vfr;

        fi_lung = (qLeft * fiLeft + qRight * fiRight) / (qRight + qLeft);
        break;
    }

}

//------------------------------------------------------------
void D0_transport::prescribe_lung_fi(solver_lumped& lum_mod){
	fi_lung = 1.;
}

//--------------------------------------------------------------
void D0_transport::initialization(){
    switch(this-> LType){
    case Perif0D:
        fi_arteriole_start.clear();
        fi_arteriole_end.clear();
        fi_capillary_start.clear();
        fi_capillary_end.clear();
        fi_venulare_start.clear();
        fi_venulare_end.clear();
        fi_vein_start.clear();
        fi_vein_end.clear();

        fi_arteriole.clear();
        fi_capillary.clear();
        fi_venulare.clear();
        fi_vein.clear();
        fi_old_arteriole.clear();
        fi_old_capillary.clear();
        fi_old_venulare.clear();
        fi_old_vein.clear();
        break;

    case Heart0D:
        fi_RA_save.clear();
        fi_RV_save.clear();
        fi_LA_save.clear();
        fi_LV_save.clear();
        fi_pul_art_start.clear();
        fi_pul_art_end.clear();
        fi_pul_vein_start.clear();
        fi_pul_vein_end.clear();
	    break;

}

}


//--------------------------------------------------------------
void D0_transport::save_variables(){

    switch (this-> LType){
    case Perif0D:
		fi_arteriole_start.push_back(fi_arteriole[0]);
        fi_arteriole_end.push_back(fi_arteriole.back());
        fi_capillary_start.push_back(fi_capillary[0]);
        fi_capillary_end.push_back(fi_capillary.back());
        fi_venulare_start.push_back(fi_venulare[0]);
        fi_venulare_end.push_back(fi_venulare.back());
        fi_vein_start.push_back(fi_vein[0]);
        fi_vein_end.push_back(fi_vein.back());
		break;

	case Heart0D:
		//nodes
		fi_RA_save.push_back(fi_RA);
		fi_RV_save.push_back(fi_RV);
		fi_LA_save.push_back(fi_LA);
		fi_LV_save.push_back(fi_LV);

		//virtual 1D
		fi_pul_art_start.push_back(fi_pul_art[0]);
		fi_pul_art_end.push_back(fi_pul_art.back());
		fi_pul_vein_start.push_back(fi_pul_vein[0]);
		fi_pul_vein_end.push_back(fi_pul_vein.back());
		break;
	}

}


//--------------------------------------------------------------
void D0_transport::save_results(string fn, const vector<double>& time, string model_name){
	string file_name;
	mkdir(("results/" + fn + "/" + model_name).c_str(),0777);

	switch (this-> LType){
    case Perif0D:
    	file_name = "results/" + fn + "/" + model_name + "/arteriole.txt";
    	//cout<<file_name;
		save_vector(file_name, fi_arteriole_start, fi_arteriole_end, time);

		file_name = "results/" + fn + "/" + model_name + "/capillary.txt";
		save_vector(file_name, fi_capillary_start, fi_capillary_end, time);

		file_name = "results/" + fn + "/" + model_name + "/venulare.txt";
		save_vector(file_name, fi_venulare_start, fi_venulare_end, time);

		file_name = "results/" + fn + "/" + model_name + "/vein.txt";
		save_vector(file_name, fi_vein_start, fi_vein_end, time);
		break;

	case Heart0D:
		file_name = "results/" + fn + "/" + model_name + "/RA.txt";
		save_vector(file_name, fi_RA_save, time);

		file_name = "results/" + fn + "/" + model_name + "/RV.txt";
		save_vector(file_name, fi_RV_save, time);

		file_name = "results/" + fn + "/" + model_name + "/LA.txt";
		save_vector(file_name, fi_LA_save, time);

		file_name = "results/" + fn + "/" + model_name + "/LV.txt";
		save_vector(file_name, fi_LV_save, time);

		file_name = "results/" + fn + "/" + model_name + "/pul_art.txt";
		save_vector(file_name, fi_pul_art_start, fi_pul_art_end, time);

		file_name = "results/" + fn + "/" + model_name + "/pul_vein.txt";
		save_vector(file_name, fi_pul_vein_start, fi_pul_vein_end, time);
		break;
	}
}


//--------------------------------------------------------------
void D0_transport::save_vector(string folder_name, const vector<double>& st, const vector<double>& en, const vector<double>& time){
    FILE *out_file = fopen(folder_name.c_str(),"w");
    //if(!out_file){cout<<"alma";}

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
    //if(!out_file){cout<<"alma";}

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
};
