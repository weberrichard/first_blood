#include "solver_lumped.h"

//--------------------------------------------------------------
solver_lumped::solver_lumped(string a_name, string a_folder)
{
	name = a_name;
	input_folder_path = a_folder;
}

solver_lumped::~solver_lumped(){}

//--------------------------------------------------------------
void solver_lumped::initialization()
{
	// setting sizes
	number_of_nodes = nodes.size();
	number_of_edges = edges.size();
	number_of_master = boundary_model_node.size();

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
	time.push_back(0.);

	// setting Eigen vars
	int nm = number_of_edges + number_of_nodes + number_of_master + 2*number_of_elastance;
	A = MatrixXd::Zero(nm,nm);
	b = VectorXd::Zero(nm);

	// building model
	build_system();
}

//--------------------------------------------------------------
vector<vector<double> > solver_lumped::solve_one_step(double dt, vector<vector<double> > boundary_coefficients)
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
	double E_act = elastance(time_act);

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
		/*int j;
		if(nodes[idx]->edge_in.size()!=0)
		{
			j = nodes[idx]->edge_in[0];
		}
		else
		{
			j = nodes[idx]->edge_out[0];
		}*/

		A(n+m+2*l+i,n+m+2*l+i) = boundary_coefficients[i][0]*1.e-6; // m3/s to ml
		A(n+m+2*l+i,m+idx) = boundary_coefficients[i][1]*mmHg_to_Pa; // Pa to mmHg
		b(n+m+2*l+i) = boundary_coefficients[i][2];

		// kinda works...
		//A(n+m+i,n+m+i) = boundary_coefficients[i][0]*1.e-6; // m3/s to ml
		//A(n+m+i,m+idx) = boundary_coefficients[i][1]*mmHg_to_Pa; // Pa to mmHg
		//b(n+m+i) = boundary_coefficients[i][2];

		// this would be in SI units
		//A(n+m+i,n+m+i) = boundary_coefficients[i][0];
		//A(n+m+i,m+idx) = boundary_coefficients[i][1];
		//b(n+m+i) = boundary_coefficients[i][2];

		// adding +1/-1 to node equation
		if(boundary_coefficients[i][1] > 0.) // ingoing edge
		{
			A(m+idx,n+m+2*l+i) = 1.;
			//A(m+idx,n+m+i) = 1.;
		}
		else // outgoing edge
		{
			A(m+idx,n+m+2*l+i) = -1.;
			//A(m+idx,n+m+i) = -1.;
		}
	}

	// actually solving the equations
	VectorXd x = A.colPivHouseholderQr().solve(b);
	//VectorXd x = A.fullPivLu().solve(b);
	//VectorXd x = A.completeOrthogonalDecomposition().solve(b);
	//VectorXd x = A.bdcSvd(ComputeThinU | ComputeThinV).solve(b);
	//VectorXd x = A.fullPivHouseholderQr().solve(b);
	//VectorXd x = A.inverse()*b;

	/*cout << A << endl;
	cout << endl << b << endl;
	cout << endl << x << endl;
	cin.get();
	cout << "rel error: " << (A*x-b).norm() / b.norm() << endl;*/

	// putting back the outputs
	for(int i=0; i<number_of_edges; i++)
	{
		edges[i]->vfr = x(i);
		edges[i]->volume_flow_rate.push_back(x(i)*1.e-6);
	}
	for(int i=0; i<number_of_nodes; i++)
	{
		nodes[i]->p = x(m+i);
		nodes[i]->y = x(m+i)/E_act;
		nodes[i]->pressure.push_back(x(m+i)*mmHg_to_Pa);
	}

	// outputs for moc model
	vector<vector<double> > out;
	for(int i=0; i<boundary_indices.size(); i++)
	{
		int idx=boundary_indices[i][3];
		double q = x(n+m+2*l+i)*1e-6;
		//double q = x(n+m+i);
		double p = x(m+idx)*mmHg_to_Pa;
		//double p = x(m+idx);
		vector<double> v{q,p};
		out.push_back(v);
	}

	return out;
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
		edges[i]->vfr_ini_non_SI = edges[i]->volume_flow_rate_initial/1.e6;
		if(edges[i]->type_code == 0) // resistance
		{
			edges[i]->par_non_SI = edges[i]->parameter/mmHg_to_Pa/1.e6;
		}
		else if(edges[i]->type_code == 1) // capacitor
		{
			edges[i]->par_non_SI = edges[i]->parameter*mmHg_to_Pa*1.e6;
		}
		else if(edges[i]->type_code == 3) // inductor
		{
			edges[i]->par_non_SI = edges[i]->parameter/mmHg_to_Pa/1.e6;
		}
		else if(edges[i]->type_code == 4) // voltage
		{
			edges[i]->par_non_SI = edges[i]->parameter/mmHg_to_Pa;
		}
		else if(edges[i]->type_code == 5) // diode
		{
			edges[i]->par_non_SI = edges[i]->parameter/mmHg_to_Pa/1.e6;
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
	// normalized version
	double tn = t * heart_rate/60.;

	// making the elastance periodic
	while(tn>1.)
	{
		tn -= 1.;
	}

	double En = 17.4073 * pow(tn,1.9) / (1.+11.2305*pow(tn,1.9)) * 1. / (1.+1.6658e7*pow(tn,21.9));

	//double En = 1.55*pow(tn/0.7,1.9)/(1.+pow(tn/0.7,1.9)) * (1./(1.+pow(tn/1.17,21.9)));

	double E = (elastance_max-elastance_min)*En + elastance_min;

	// E = E*mmHg_to_Pa*1.e6; // mmHg/ml to SI: Pa/m3 

	return E;
}

//--------------------------------------------------------------
double solver_lumped::elastance_derived(double t)
{
	// normalized version
	double tn = t * heart_rate/60.;

	// making the elastance periodic
	while(tn>1.)
	{
		tn -= 1.;
	}

	double Enp = (9.450202509727443e-16*pow(tn,0.9) - 1.6570681411267346e-7*pow(tn,22.8) - 2.037762561602155e-6*pow(tn,24.7))/(pow(0.0890432 + pow(tn,1.9),2.)*pow(6.003121623244087e-8 + pow(tn,21.9),2.));

	double Ep = (elastance_max-elastance_min)*Enp;

	//Ep = Ep*mmHg_to_Pa*1.e6;

	return Ep;
}
