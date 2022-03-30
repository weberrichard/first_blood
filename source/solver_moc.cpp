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

	// counting the number of division points
	sum_division_points = 0;
	for(unsigned int i=0; i<number_of_edges; i++)
	{
		sum_division_points += edges[i]->division_points;
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
double solver_moc::timesteps(int &idx)
{
	// getting the real time step i.e. the minimum timestep from every edge
	for(unsigned int j=0; j<forward_edges.size(); j++)
	{
		edges[forward_edges[j]]->new_timestep();
	}

	// getting the min timestep
	double t_min = min_time(idx);

	return t_min;
}

//--------------------------------------------------------------
void solver_moc::boundaries(int e_idx, double t_act)
{
	// find all neigbouring nodes
	vector<int> node_idx{edges[e_idx]->node_index_start,edges[e_idx]->node_index_end};

	for(unsigned int i=0; i<node_idx.size(); i++)
	{
		if(nodes[node_idx[i]]->is_master_node == false)
		{
			if(nodes[node_idx[i]]->is_upstream_boundary) // handling the upstream boundary
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

					double t = t_act-period*(time_upstream.back()-time_upstream[0]);

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
				double v_h = value_upstream[index_upstream+1]; // pressure at higher index
				double v_l = value_upstream[index_upstream]; // pressure at lower index
				double t_h = time_upstream[index_upstream+1]; // time at higher index
				double t_l = time_upstream[index_upstream]; // time at lower index

				double t_in = t_act-period*time_upstream.back(); // actual time of the simulation
				double v_in = (v_h-v_l)/(t_h-t_l) * (t_in-t_l) + v_l; // actual pressure of the simulation

				double q_in=0., p_in;
				// there might be several outgoing edge from an upstream node
				for(unsigned int j=0; j<nodes[node_idx[i]]->edge_out.size(); j++)
				{
					int ei = nodes[node_idx[i]]->edge_out[j];
					// calculating vp and pp values of edge
					if(type_upstream == 0) // pressure BC
					{
						double dt = edges[ei]->dt_act;
						p_in = v_in; // v_in is pressure
						q_in += edges[ei]->upstream_boundary_p(dt, p_in);
					}
					else if(type_upstream == 1) // vfr BC
					{
						double dt = edges[ei]->dt_act;
						p_in = edges[ei]->upstream_boundary_q(dt, v_in);
						q_in += v_in; // v_in is volume flow rate
					}
					else if(type_upstream == 2) // v BC
					{
						double dt = edges[ei]->dt_act;
						double q=0.;
						p_in = edges[ei]->upstream_boundary_v(dt, v_in, q);
						q_in += q; // v_in is volume flow rate
					}
				}

				if(nodes[node_idx[i]]->do_save_memory)
				{
					nodes[node_idx[i]]->save_field_variables(t_act,p_in,q_in);
				}
			}
			else if(nodes[node_idx[i]]->type_code == 1) // perifera
			{
				// there can be only one incoming edge
				int edge_index = nodes[node_idx[i]]->edge_in[0];
				double dt = edges[edge_index]->dt_act;
				double p_out = nodes[node_idx[i]]->pressure_out;
				double q_out = edges[edge_index]->boundary_periferia(dt,p_out);

				if(nodes[node_idx[i]]->do_save_memory)
				{
					nodes[node_idx[i]]->save_field_variables(t_act,p_out,q_out);
				}
			}
			else if(nodes[node_idx[i]]->type_code == 0) // junctions
			{
				// first the node pressure is calculated

				// temp variables
				double num=0., denum=0.;

				// leakage part if there is, otherwise nb=[0,0]
				vector<double> nb = nodes[node_idx[i]]->boundary_coefficients();
				num += nb[0];
				denum += nb[1];

				// for saving edge coefs
				vector<double> edge_coefs_in(2,0.), edge_coefs_out(2,0.);
				// ingoing edges
				for(unsigned int j=0; j<nodes[node_idx[i]]->edge_in.size(); j++)
				{
					int edge_index = nodes[node_idx[i]]->edge_in[j];
					double dt = t_act - edges[edge_index]->time.back();
					vector<double> eb = edges[edge_index]->boundary_end_coefficients(dt);
					num   -= eb[1];
					denum += eb[0];
					if(edge_index == e_idx)
					{
						edge_coefs_in = eb;
					}
				}

				// outgoing edges
				for(unsigned int j=0; j<nodes[node_idx[i]]->edge_out.size(); j++)
				{
					int edge_index = nodes[node_idx[i]]->edge_out[j];
					double dt = t_act - edges[edge_index]->time.back();
					vector<double> eb = edges[edge_index]->boundary_start_coefficients(dt);
					num   += eb[1];
					denum -= eb[0];
					if(edge_index == e_idx)
					{
						edge_coefs_out = eb;
					}
				}

				// nodal pressure
				double p_nodal = num/denum;

				nodes[node_idx[i]]->boundary_variables(p_nodal,t_act);

				// edge velocity and pressure
				if(i==1) // end node
				{
					double q = edge_coefs_in[0]*p_nodal + edge_coefs_in[1];
					//double dt = t_act - edges[e_idx]->time.back();
					double dt = edges[e_idx]->dt_act;
					edges[e_idx]->boundary_end_variables(dt,p_nodal,q);
				}
				else // start node
				{
					double q = edge_coefs_out[0]*p_nodal + edge_coefs_out[1];
					//double dt = t_act - edges[e_idx]->time.back();
					double dt = edges[e_idx]->dt_act;
					edges[e_idx]->boundary_start_variables(dt,p_nodal,q);
				}
			}
		}
	}
}

//--------------------------------------------------------------
double solver_moc::min_time(int &idx)
{
	double t_min = 1.e10;
	idx = -1;
	for(unsigned int j=0; j<forward_edges.size(); j++)
	{
		if(t_min > edges[forward_edges[j]]->time.back()+edges[forward_edges[j]]->dt_act)
		{
			t_min = edges[forward_edges[j]]->time.back()+edges[forward_edges[j]]->dt_act;
			idx = forward_edges[j];
		}
	}

	return t_min;
}

//--------------------------------------------------------------
void solver_moc::post_process()
{
	for(unsigned int j=0; j<forward_edges.size(); j++)
	{
		// interpolating to equidistant mesh
		edges[forward_edges[j]]->interpolate();
		// updating every field variables
		edges[forward_edges[j]]->update_variables();
		// saveing start and end points vars in time
		if(edges[forward_edges[j]]->do_save_memory)
		{
			edges[forward_edges[j]]->save_field_variables();
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
void solver_moc::convert_time_series()
{
	if(type_upstream == 0)
	{
		for(unsigned int i=0; i<value_upstream.size(); i++)
		{
			value_upstream[i] *= mmHg_to_Pa;
			value_upstream[i] += atmospheric_pressure;
		}
	}
	else if(type_upstream == 1)
	{
		for(unsigned int i=0; i<value_upstream.size(); i++)
		{
			value_upstream[i] *= 1.e-6;
		}
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
		if(edge_id.compare(edges[i]->ID) == 0)
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

//--------------------------------------------------------------
vector<int> solver_moc::edge_to_node(vector<int> edge_idx)
{
	vector<bool> node_bool(number_of_nodes,false);
	for(unsigned int i=0; i<edge_idx.size(); i++)
	{
		node_bool[edges[i]->node_index_start] = true;
		node_bool[edges[i]->node_index_end] = true;
	}

	vector<int> node_idx;
	for(unsigned int i=0; i<number_of_nodes; i++)
	{
		if(node_bool[i])
		{
			node_idx.push_back(i);
		}
	}
	return node_idx;
}

//--------------------------------------------------------------
void solver_moc::clear_save_memory()
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
void solver_moc::set_save_memory(vector<string> edge_list, vector<string> node_list)
{
	for(int i=0; i<edge_list.size(); i++)
	{
		int idx = edge_id_to_index(edge_list[i]);
		edges[idx]->do_save_memory = true;
	}
	for(int i=0; i<node_list.size(); i++)
	{
		int idx = node_id_to_index(node_list[i]);
		nodes[idx]->do_save_memory = true;
	}
}

