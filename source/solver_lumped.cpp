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
}

//--------------------------------------------------------------
void solver_lumped::set_newton_size()
{
	number_of_moc = boundary_indices.size();

	// setting Eigen vars for nonlinear solvr
	int N = number_of_edges + number_of_nodes + 2*number_of_elastance + 2*number_of_moc;
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
void solver_lumped::initialization_newton()
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

}

//--------------------------------------------------------------
/*vector<vector<double> > solver_lumped::solve_one_step(double dt, vector<vector<double> > boundary_coefficients)
{
	// increasing time
	double time_prev = time.back();
	double time_act = time_prev+dt;
	time.push_back(time_act);

	// sizes of nodes and edges
	int n=number_of_nodes, m=number_of_edges, l=number_of_elastance;

	// tracing the virtual nodes of elastance
	int i_elas=0;

	// actual elastance
	double E_act = elastance(time_act,edges[i]->parameter);

	// edges
	for(int i=0; i<number_of_edges; i++)
	{
		int i1 = edges[i]->node_index_start;
		int i2 = edges[i]->node_index_end;

		if(edges[i]->type_code == 0) // resistor
		{
			A(i,m+i2) = 1.;
			A(i,m+i1) = -1.;
			A(i,i) = edges[i]->par_non_SI; // R
			//A(i,i) = edges[i]->parameter; // R
		}
		else if(edges[i]->type_code == 1) // capacitor
		{
			A(i,m+i2) = 1.;
			A(i,m+i1) = -1.;
			A(i,i) = dt/edges[i]->par_non_SI; // dt/C
			b(i) = nodes[i2]->p - nodes[i1]->p; // pj-pi from previous timestamp
		}
		else if(edges[i]->type_code == 2) // elastance
		{
			// basic equation for the edge
			A(i,m+n+i_elas+1) = 1.;
			A(i,m+n+i_elas) = -1.;
			A(i,i) = dt;
			b(i) = nodes[i2]->y-nodes[i1]->y;

			// equations for the virtual nodes
			A(n+m+i_elas,m+i1) = -1.;
			A(n+m+i_elas+1,m+i2) = -1.;
			A(n+m+i_elas,n+m+i_elas) = E_act;
			A(n+m+i_elas+1,n+m+i_elas+1) = E_act;
			i_elas+=2;
		}
		else if(edges[i]->type_code == 3) // inductor
		{
			A(i,m+i2) = 1.;
			A(i,m+i1) = -1.;
			A(i,i) = edges[i]->par_non_SI/dt; // L/dt
			b(i) = edges[i]->par_non_SI/dt * edges[i]->vfr; // L/dt * qk
		}
		else if(edges[i]->type_code == 4) // voltage source
		{
			A(i,m+i2) = 1.;
			A(i,m+i1) = -1.;
			b(i) = edges[i]->par_non_SI; // p_const
		}
		else if(edges[i]->type_code == 5) // diode
		{
			double p2 = nodes[i2]->p;
			double p1 = nodes[i1]->p;
			bool is_closed;

			if(p1>p2) // diode is open
			{
				A(i,m+i2) = 1.;
				A(i,m+i1) = -1.;
				A(i,i) = edges[i]->par_non_SI; // R
				is_closed = false;
			}
			else // diode is closed
			{
				A(i,m+i2) = 1.;
				A(i,m+i1) = -1.;
				A(i,i) = 1.e10*edges[i]->par_non_SI; // R*1.e10
				is_closed = true;
			}
		}
	}

	// nodes
	for(int i=0; i<number_of_nodes; i++)
	{
		if(nodes[i]->is_ground == false) // non-ground nodes, intersections
		{
			for(int j=0; j<nodes[i]->edge_in.size(); j++)
			{
				A(m+i,nodes[i]->edge_in[j]) = 1;
			}
			for(int j=0; j<nodes[i]->edge_out.size(); j++)
			{
				A(m+i,nodes[i]->edge_out[j]) = -1;
			}
		}
		else // ground nodes, pi = p0[mmHg]
		{
			A(m+i,m+i) = 1;
			b(m+i) = 1.e5/mmHg_to_Pa;
		}
	}

	// master boundary node from MOC
	for(int i=0; i<boundary_indices.size(); i++)
	{
		// setting indecies
		int idx = boundary_indices[i][3];

		A(n+m+2*l+i,n+m+2*l+i) = boundary_coefficients[i][0]*1.e-6; // m3/s to ml
		A(n+m+2*l+i,m+idx) = boundary_coefficients[i][1]*mmHg_to_Pa; // Pa to mmHg
		b(n+m+2*l+i) = boundary_coefficients[i][2];

		// adding +1/-1 to node equation
		if(boundary_coefficients[i][1] > 0.) // ingoing edge
		{
			A(m+idx,n+m+2*l+i) = 1.;
		}
		else // outgoing edge
		{
			A(m+idx,n+m+2*l+i) = -1.;
		}
	}

	// actually solving the equations
	VectorXd x = A.colPivHouseholderQr().solve(b);

	// putting back the outputs
	for(int i=0; i<number_of_edges; i++)
	{
		edges[i]->vfr = x(i);
		if(edges[i]->do_save_memory)
		{
			edges[i]->volume_flow_rate.push_back(x(i)*1.e-6);
		}
	}
	for(int i=0; i<number_of_nodes; i++)
	{
		nodes[i]->p = x(m+i);
		nodes[i]->y = x(m+i)/E_act;
		if(nodes[i]->do_save_memory)
		{
			nodes[i]->pressure.push_back(x(m+i)*mmHg_to_Pa);
		}
	}

	// outputs for moc model
	vector<vector<double> > out;
	for(int i=0; i<boundary_indices.size(); i++)
	{
		int idx=boundary_indices[i][3];
		double q = x(n+m+2*l+i)*1e-6;
		double p = x(m+idx)*mmHg_to_Pa;
		vector<double> v{q,p};
		out.push_back(v);
	}

	return out;
}
*/

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

