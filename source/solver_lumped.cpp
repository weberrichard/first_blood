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
				Jac(i,i) = par; // R
				f(i) = x(m+i2) - x(m+i1) + par*x(i);
			}
			else // diode is closed
			{
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
int NX(double L,double dx, int N) {
    if (floor(L / dx) > N) {
        return N;
    }
    else {
        return floor(L / dx);
    }
}


//-----------------------------------------
void Virt1DforLum(vector<double> &fi_old, vector<double> &fi, double v, double dt, double dx, int n, double fiStartNode, double fiEndNode) {
    vector<double> fi_tmp = fi;

    for (int i = 1; i < n - 1; i++) {

        if (v > 0) {
            fi[i] = fi_old[i] - v * dt / dx * (fi_old[i] - fi_old[i - 1]);
        }
        else {
            fi[i] = fi_old[i] - v * dt / dx * (fi_old[i + 1] - fi_old[i]);
        }
    }

    //BC
    if (v > 0) {

        fi[n - 1] = fi_old[n - 1] - v * dt / dx * (fi_old[n - 1] - fi_old[n - 2]);
        fi[0] = fiStartNode;
    }
    else {
        fi[n - 1] = fiEndNode; 
        fi[0] = fi_old[0] - v * dt / dx * (fi_old[1] - fi_old[0]);
    }

    fi_old = fi_tmp;

}


//--------------------------------
D0Transport::D0Transport(LumpedType LType, vector<string> sv, TransportType TType):LType(LType), TType(TType) {
    switch (LType) {
    case PerifCoronary0D:
    //no idea
        break;
    case Perif0D:

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
        fi_capillary.assign(nx_arteriole, 0.);
        fi_old_capillary.assign(nx_arteriole, 0.);

        //venulare
        nx_venulare = NX(L_venulare, stod(sv[11],0), 5);
        dx_venulare = L_venulare / (nx_venulare - 1);
        fi_venulare.assign(nx_arteriole, 0.);
        fi_old_venulare.assign(nx_arteriole, 0.);

        //vein
        nx_vein = NX(L_vein, stod(sv[12],0), 30);
        dx_vein = L_vein / (nx_vein - 1);
        fi_vein.assign(nx_arteriole, 0.);
        fi_old_vein.assign(nx_arteriole, 0.);


        break;
    case Heart0D:


        break;
    default:
        cout << "Unknown 0D type for transport.";
    }
}


//-----------------------------------
void D0Transport::UpdateFi(int LumIndex, double dt, double masterFi, vector<solver_lumped*> lum) {
	//index of the heart model
	int heartIndex = -1;// this could be done only at the initialization i guess...
        for (int i=0;i<lum.size(); i++){
        	if (lum[i]->name == "heart_kim_lit"){ heartIndex = i;}
        }
        if(heartIndex == -1){ cout<<"???"; return; }


    switch (this-> LType) {
    case PerifCoronary0D:
        //no idea
        break;
    case Perif0D:
      
        Virt1DforLum(fi_old_arteriole, fi_arteriole, lum[LumIndex]->edges[1]->vfr * ml_to_m3 / A_arteriole, dt, dx_arteriole, nx_arteriole, lum[LumIndex]->nodes[1]->RBC_fi0Dn, lum[LumIndex]->nodes[2]->RBC_fi0Dn);
        Virt1DforLum(fi_old_capillary, fi_capillary, lum[LumIndex]->edges[2]->vfr * ml_to_m3 / A_capillary, dt, dx_capillary, nx_capillary, lum[LumIndex]->nodes[2]->RBC_fi0Dn, lum[LumIndex]->nodes[3]->RBC_fi0Dn);
        Virt1DforLum(fi_old_venulare, fi_venulare, lum[LumIndex]->edges[3]->vfr * ml_to_m3 / A_venulare, dt, dx_venulare, nx_venulare, lum[LumIndex]->nodes[3]->RBC_fi0Dn, lum[LumIndex]->nodes[4]->RBC_fi0Dn);
        Virt1DforLum(fi_old_vein, fi_vein, lum[LumIndex]->edges[4]->vfr * ml_to_m3 / A_vein, dt, dx_vein, nx_vein, lum[LumIndex]->nodes[4]->RBC_fi0Dn, lum[heartIndex]->nodes[0]->RBC_fi0Dn); // n1 is the master node
            
        //nodes
        // not sure if the virtual 1D or the nodes should be updated first
        UpdatePerifLumNodes(1, LumIndex, masterFi , fi_arteriole[0], lum);// masterFi is the master node's fi
        UpdatePerifLumNodes(2, LumIndex, fi_arteriole.back(), fi_capillary[0], lum);
        UpdatePerifLumNodes(3, LumIndex, fi_capillary.back(), fi_venulare[0], lum);
        UpdatePerifLumNodes(4, LumIndex, fi_venulare.back(), fi_vein[0], lum);

        //master node
        //node connecting heart and perifs


        break;
    case Heart0D:// that will be fun...
        //right atrium


        break;
    default:
        cout << "Unknown 0D type for transport.";
    }
    return;
}

//------------------------------------
//updates fi parameters of perif nodes
void D0Transport::UpdatePerifLumNodes(int LumNodeIndex, int LumIndex, double fiLeft, double fiRight, vector<solver_lumped*> lum) {
    int a, b, c;

    //1 if q flows towards the node
    if (lum[LumIndex]->edges[LumNodeIndex]->vfr < 0) { a = 1; } //Resistance edge
    if (lum[LumIndex]->edges[LumNodeIndex - 1]->vfr > 0) { b = 1; } //Inductance edge
    if (lum[LumIndex]->edges[LumNodeIndex + 4]->vfr < 0) { c = 1; } //Capacitance edge

    double qLeft, qRight, qDown;

    switch (a*4 + b*2 + c) {
    case 1:
        lum[LumIndex]->nodes[LumNodeIndex]->RBC_fi0Dn = fiRight;
        break;

    case 2:
        lum[LumIndex]->nodes[LumNodeIndex]->RBC_fi0Dn = fiLeft;
        break;

    case 3:
        qLeft = lum[LumIndex]->edges[LumNodeIndex - 1]->vfr;
        qDown = lum[LumIndex]->edges[LumNodeIndex - 4]->vfr;
        lum[LumIndex]->nodes[LumNodeIndex]->RBC_fi0Dn = (qLeft*fiLeft + qDown *fiRight)/(qDown + qLeft);
        break;

    case 4:
        lum[LumIndex]->nodes[LumNodeIndex]->RBC_fi0Dn = fiRight; //same az case 1
        break;

    case 5:
        lum[LumIndex]->nodes[LumNodeIndex]->RBC_fi0Dn = fiRight; //both has the same fi
        break;

    case 6:
        qLeft = lum[LumIndex]->edges[LumNodeIndex - 1]->vfr; // in ml/s, the dimension does not matter
        qRight = lum[LumIndex]->edges[LumNodeIndex]->vfr;
        lum[LumIndex]->nodes[LumNodeIndex]->RBC_fi0Dn = (qLeft * fiLeft + qRight * fiRight) / (qRight + qLeft);
        break;
    }
}


