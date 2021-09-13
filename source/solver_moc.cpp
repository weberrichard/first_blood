#include "solver_moc.h"

//--------------------------------------------------------------
solver_moc::solver_moc(string a_name, string a_folder)
{
	name = a_name;
	input_folder_path = a_folder;
}

solver_moc::~solver_moc(){}

//--------------------------------------------------------------
void solver_moc::initialization(double p_init)
{
	// setting the size of nodes and edges
	number_of_nodes = nodes.size();
	number_of_edges = edges.size();

	for(unsigned int i=0; i<number_of_nodes; i++)
	{
		nodes[i]->initialization(p_init);
	}
	for(unsigned int i=0; i<number_of_edges; i++)
	{
		edges[i]->initialization(p_init);
	}

	// setting back the pressure_upstream interpolation index to 0
	index_upstream = 0;

	// setting back the period counter
	period = 0;

	// also storing the initial pressure
	pressure_initial = p_init;

	// setting back time
	time.clear();

	if(time_upstream.size()>0)
	{
		time.push_back(time_upstream[0]);
	}
	else
	{
		time.push_back(0.);
	}

	// building up the topology
	build_system();
}

//--------------------------------------------------------------
void solver_moc::build_system()
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
void solver_moc::full_solver(string node_id, double time_end)
{
	forward_solver(node_id, time_end);

	int node_idx = node_id_to_index(node_id);
	while(nodes[node_idx]->type_code != 2) // checking if it is upstream boundary
	{
		// performing backward calculation
		int edge_idx = backward_tree(node_id);
		// boundary conditions from last forward calculation
		vector<double> t_d = time;
		vector<double> p_d = nodes[node_idx]->pressure;
		vector<double> vfr_d(t_d.size(),0.);
		for(unsigned int i=0; i<vfr_d.size(); i++)
		{  
			for(unsigned int j=0; j<nodes[node_idx]->edge_out.size(); j++)
			{
				vfr_d[i] += edges[nodes[node_idx]->edge_out[j]]->volume_flow_rate_start[i];
			}
		}

		vector<vector<double> > A_u = edges[edge_idx]->backward_solver(t_d, p_d, vfr_d);

		// new boundary conditions
		time_upstream = A_u[0];
		pressure_upstream = A_u[1];
		time_end = time_upstream.back();

		// new upstream node
		node_id = edges[edge_idx]->node_name_start;
		node_idx = edges[edge_idx]->node_index_start;

		// forward calculation
		forward_solver(node_id, time_end);
	}
}

//--------------------------------------------------------------
void solver_moc::forward_solver(string node_id, double time_end)
{
	// setting initial conditions
	initialization(pressure_initial);

	// getting the forward tree edges and nodes
	forward_tree(node_id);

	// starting the main solver cycle
	int i=0;
	while(time.back() < time_end)
	{
		// updating on timestep
		double dt = solve_one_step();

		// interpolation, updating and saving field variables
		post_process(dt);

		i++;
	}

	// saving the final number of time steps
	number_of_timesteps = time.size();
}

//--------------------------------------------------------------
double solver_moc::solve_one_step()
{
	// getting the real time step i.e. the minimum timestep from every edge
	double dt_real=1.e10;
	for(unsigned int j=0; j<forward_edges.size(); j++)
	{
		double dt=edges[forward_edges[j]]->new_timestep();
		if(dt<dt_real)
		{
			dt_real = dt;
		}
	}
	if(dt_real<0.)
	{
		printf("\n !WARNING! time step is negative: %6.3e during FORWARD calculation", dt_real);
		cout << endl;
	}

	double new_time = time.back() + dt_real;
	//if(new_time > time_end)
	//{
	//	dt_real = time_end - time.back();
	//	new_time = time_end;
	//}
	time.push_back(new_time);

	// calculating new pressure and velocity field in inner points
	for(unsigned int j=0; j<forward_edges.size(); j++)
	{
		edges[forward_edges[j]]->solve();
	}

	// handling nodes i.e. boundaries conditions of edges
	boundaries(dt_real);

	return dt_real;
}

//--------------------------------------------------------------
void solver_moc::post_process(double dt)
{
	for(unsigned int j=0; j<forward_edges.size(); j++)
	{
		// interpolating to equidistant mesh
		edges[forward_edges[j]]->interpolate(dt);
		// updating every field variables
		edges[forward_edges[j]]->update_variables(dt);
		// saveing start and end points vars in time
		edges[forward_edges[j]]->save_field_variables();
	}
}

//--------------------------------------------------------------
void solver_moc::boundaries(double dt)
{
	for(unsigned int i=0; i<forward_nodes.size(); i++)
	{
		if(nodes[forward_nodes[i]]->is_master_node == false)
		{
			if(nodes[forward_nodes[i]]->is_upstream_boundary) // handling the upstream boundary
			{
				// finding the position for linear interpolation
				int j=index_upstream;
				bool got_it = false;
				while(!got_it)
				{
					// making the inlet function periodic
					if(j > time_upstream.size()-1)
					{
						j -= time_upstream.size();
						period += 1;
					}

					double t = time.back()-period*(time_upstream.back()-time_upstream[0]);

					if(t >= time_upstream[j] && t <= time_upstream[j+1])
					{
						got_it=true;
						index_upstream = j;
					}
					else
					{
						j++;
					}
				}

				// interpolating
				double p_h = pressure_upstream[index_upstream+1]; // pressure at higher index
				double p_l = pressure_upstream[index_upstream]; // pressure at lower index
				double t_h = time_upstream[index_upstream+1]; // time at higher index
				double t_l = time_upstream[index_upstream]; // time at lower index

				double t_in = time.back()-period*time_upstream.back(); // actual time of the simulation
				double p_in = (p_h-p_l)/(t_h-t_l) * (t_in-t_l) + p_l; // actual pressure of the simulation

				double q_in=0.;
				// there might be several outgoing edge from an upstream node
				for(unsigned int j=0; j<nodes[forward_nodes[i]]->edge_out.size(); j++)
				{
					int edge_index = nodes[forward_nodes[i]]->edge_out[j];
					// calculating vp and pp values of edge
					q_in += edges[edge_index]->upstream_boundary(dt, p_in);
				}

				nodes[forward_nodes[i]]->pressure.push_back(p_in);
				nodes[forward_nodes[i]]->volume_flow_rate.push_back(q_in);
			}
			else if(nodes[forward_nodes[i]]->type_code == 1) // perifera
			{
				// there can be only one incoming edge
				int edge_index = nodes[forward_nodes[i]]->edge_in[0];
				double p_out = nodes[forward_nodes[i]]->pressure_out;
				double q_out = edges[edge_index]->boundary_periferia(dt,p_out);

				nodes[forward_nodes[i]]->pressure.push_back(p_out);
				nodes[forward_nodes[i]]->volume_flow_rate.push_back(q_out);
			}
			else if(nodes[forward_nodes[i]]->type_code == 0) // junctions
			{
				// first the node pressure is calculated

				// temp variables
				double num=0., denum=0.;

				// leakage part if there is, otherwise nb=[0,0]
				vector<double> nb = nodes[forward_nodes[i]]->boundary_coefficients();
				num += nb[0];
				denum += nb[1];

				// handling edges

				// for saving edge coefs
				vector<vector<double> > edge_coefs_in, edge_coefs_out;

				// ingoing edges
				for(unsigned int j=0; j<nodes[forward_nodes[i]]->edge_in.size(); j++)
				{
					int edge_index = nodes[forward_nodes[i]]->edge_in[j];
					vector<double> eb = edges[edge_index]->boundary_end_coefficients(dt);
					num   -= eb[1];
					denum += eb[0];
					edge_coefs_in.push_back(eb);
				}

				// outgoing edges
				for(unsigned int j=0; j<nodes[forward_nodes[i]]->edge_out.size(); j++)
				{
					int edge_index = nodes[forward_nodes[i]]->edge_out[j];
					vector<double> eb = edges[edge_index]->boundary_start_coefficients(dt);
					num   += eb[1];
					denum -= eb[0];
					edge_coefs_out.push_back(eb);
				}

				// nodal pressure
				double p_nodal = num/denum;
				nodes[forward_nodes[i]]->boundary_variables(p_nodal);

				double q_nodal = 0.0;
				// edge velocity and pressure
				for(unsigned int j=0; j<nodes[forward_nodes[i]]->edge_in.size(); j++)
				{
					double q = edge_coefs_in[j][0]*p_nodal + edge_coefs_in[j][1];
					q_nodal += q; 
					int edge_index = nodes[forward_nodes[i]]->edge_in[j];
					edges[edge_index]->boundary_end_variables(dt,p_nodal,q);
				}

				for(unsigned int j=0; j<nodes[forward_nodes[i]]->edge_out.size(); j++)
				{
					double q = edge_coefs_out[j][0]*p_nodal + edge_coefs_out[j][1];
					q_nodal -= q;
					int edge_index = nodes[forward_nodes[i]]->edge_out[j];
					edges[edge_index]->boundary_start_variables(dt,p_nodal,q);
				}
			}
		}
	}
}

//--------------------------------------------------------------
void solver_moc::full_tree()
{
	forward_nodes.clear();
	forward_edges.clear();
	for(int i=0; i<number_of_nodes; i++)
	{
		forward_nodes.push_back(i);
		if(nodes[i]->is_master_node == false) // not master node
		{
			if(nodes[i]->type_code == 2) // heart
			{
				nodes[i]->is_upstream_boundary = true;
			}
		}
	}
	for(int i=0; i<number_of_edges; i++)
	{
		forward_edges.push_back(i);
	}
}

//--------------------------------------------------------------
void solver_moc::forward_tree(string node_id)
{
	// initilizing the control variables
	forward_edges.clear();
	forward_nodes.clear();
	for(unsigned int i=0; i<edges.size(); i++)
	{
		edges[i]->do_solve = false;
	}

	// finding node index
	int node_index = node_id_to_index(node_id);
	if(node_index == -1)
	{
		cout << "\n ! ERROR ! Node ID not found: " << node_id << "\n Exiting..." << endl;
		exit(-1);
	}
	// setting upstream node
	nodes[node_index]->is_upstream_boundary = true;

	// finding edges directly connected to the node
	for(unsigned int i=0; i<number_of_edges; i++)
	{
		if(node_index == edges[i]->node_index_start)
		{
			forward_edges.push_back(i);
			forward_nodes.push_back(edges[i]->node_index_start);
			forward_nodes.push_back(edges[i]->node_index_end);
		}
	}

	vector<int> ns;
	vector<int> es;
	for(unsigned int i=0; i<number_of_edges; i++)
	{
		ns.push_back(edges[i]->node_index_start);
		es.push_back(i);
	}

	// finding the rest of the forward tree
	int i=0;
	while(i<forward_edges.size())
	{
		int j=0;
		bool no_new = true;
		while(j<ns.size() && no_new)
		{
			if(edges[forward_edges[i]]->node_index_end == ns[j])
			{
				forward_edges.push_back(es[j]);
				forward_nodes.push_back(edges[es[j]]->node_index_start);
				forward_nodes.push_back(edges[es[j]]->node_index_end);
				no_new = false;
				i = 0;

				ns.erase(ns.begin()+j);
				es.erase(es.begin()+j);
			}
			j++;
		}

		if(no_new)
		{
			i++;
		}
	}

	forward_nodes = unique(forward_nodes);

	for(unsigned int i=0; i<forward_edges.size(); i++)
	{
		edges[forward_edges[i]]->do_solve = true;
	}
}

//-------------------------------------------------------------- 
vector<int> solver_moc::unique(vector<int> x)
{
	vector<int> out;
	for(int i=0; i<x.size(); i++)
	{
		bool unique = true;
		for(int j=0; j<out.size(); j++)
		{
			if(x[i] == out[j])
			{
			  unique = false;
			  break;
			}
		}
		if(unique)
		{
			out.push_back(x[i]);
		}
	}
  return out;
}

//--------------------------------------------------------------
int solver_moc::backward_tree(string node_id)
{
	int edge_index = -1;
	int i=0;
	while(edge_index == -1 && i<number_of_edges)
	{
		if(edges[i]->node_name_end == node_id)
		{
			edge_index = i;
		}
		i++;
	}

	return edge_index;
}

//--------------------------------------------------------------
void solver_moc::set_constants(double g, double rho, double nu, double mmHg, double p0, double beta)
{
	gravity = g;
	density = rho;
	kinematic_viscosity = nu;
	mmHg_to_Pa = mmHg;
	atmospheric_pressure = p0;
	for(unsigned int i=0; i<edges.size(); i++)
	{
		edges[i]->gravity = g;
		edges[i]->density = rho;
		edges[i]->kinematic_viscosity = nu;
		edges[i]->atmospheric_pressure = p0;
		edges[i]->beta = beta;
	}
	for(unsigned int i=0; i<nodes.size(); i++)
	{
		nodes[i]->density = rho;
		nodes[i]->pressure_out = atmospheric_pressure;
	}
}

//--------------------------------------------------------------
void solver_moc::convert_pt_series()
{
	for(unsigned int i=0; i<pressure_upstream.size(); i++)
	{
		pressure_upstream[i] *= mmHg_to_Pa;
		pressure_upstream[i] += atmospheric_pressure;
	}
}

//--------------------------------------------------------------
int solver_moc::node_id_to_index(string node_id)
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
		cout << "\n !!!WARNING!!!\n solver_moc::node_id_to_index function\nNode is not existing, node_id: " << node_id << "\n Continouing..." << endl;
	}
	return idx;
}

//--------------------------------------------------------------
int solver_moc::edge_id_to_index(string edge_id)
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
		cout << "\n!!!WARNING!!!\n solver_moc::edge_id_to_index function\n Node is not existing, edge_id: " << edge_id << "\n Continouing..." << endl;
	}
	return idx;
}