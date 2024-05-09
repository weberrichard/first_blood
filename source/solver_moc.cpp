#include "solver_moc.h"

//--------------------------------------------------------------
solver_moc::solver_moc(string a_name, string a_folder)
{
	name = a_name;
	input_folder_path = a_folder;
}

solver_moc::~solver_moc(){}

//--------------------------------------------------------------
void solver_moc::initialization(double p_init, int material_type, double RBC_init)
{
	// setting the size of nodes and edges
	number_of_nodes = nodes.size();
	number_of_edges = edges.size();

	for(unsigned int i=0; i<number_of_nodes; i++)
	{
		nodes[i]->initialization(p_init, RBC_init);
	}
	for(unsigned int i=0; i<number_of_edges; i++)
	{
		edges[i]->initialization(p_init,material_type, RBC_init);
	}

	// setting back the pressure_upstream interpolation index to 0
	index_upstream.assign(type_upstream.size(),0);

	// setting back the period counter
	period.assign(type_upstream.size(),0.);

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

	// setting newton iteration variables
	for(unsigned int i=0; i<number_of_edges; i++)
	{
		int n1 = nodes[edges[i]->node_index_start]->edge_in.size() + nodes[edges[i]->node_index_start]->edge_out.size();
		int n2 = nodes[edges[i]->node_index_end]->edge_in.size() + nodes[edges[i]->node_index_end]->edge_out.size();
		edges[i]->set_newton_size(n1,n2);
	}
}

//--------------------------------------------------------------
void solver_moc::initialization_newton(VectorXd &x, int N, int moc_edge_index, int edge_end)
{
	if(edge_end==1)
	{
		x(N) = edges[moc_edge_index]->initialization_newton_end();
	}
	else
	{
		x(N) = edges[moc_edge_index]->initialization_newton_start();
	}
}

//--------------------------------------------------------------
void solver_moc::substitute_newton(int moc_edge_index, int edge_end, double t_act, double p, double q)
{
	if(edge_end==1)
	{
		edges[moc_edge_index]->boundary_substitute_end(t_act, p, q);
		int node_index = edges[moc_edge_index]->node_index_end;
		nodes[node_index]->boundary_variables(p,t_act);
	}
	else
	{
		edges[moc_edge_index]->boundary_substitute_start(t_act, p, q);
		int node_index = edges[moc_edge_index]->node_index_start;
		nodes[node_index]->boundary_variables(p,t_act);
	}
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
void solver_moc::timesteps()
{
	// getting the real time step i.e. the minimum timestep from every edge
	for(unsigned int j=0; j<number_of_edges; j++)
	{
		edges[j]->new_timestep();
	}
}

//--------------------------------------------------------------
double solver_moc::min_time(int &idx)
{
	double t_min = 1.e10;
	idx = -1;
	for(unsigned int j=0; j<number_of_edges; j++)
	{
		if(t_min > edges[j]->time.back()+edges[j]->dt_act)
		{
			t_min = edges[j]->time.back()+edges[j]->dt_act;
			idx = j;
		}
	}

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
			if(nodes[node_idx[i]]->upstream_boundary>-1) // handling the upstream boundary
			{
				int up_idx = nodes[node_idx[i]]->upstream_boundary;
				// finding the position for linear interpolation
				int j=index_upstream[up_idx];
				bool got_it = false;
				while(!got_it)
				{
					// making the inlet function periodic
					if(j >= time_upstream[up_idx].size()-1)
					{
						j -= time_upstream[up_idx].size();
						period[up_idx] += 1;
					}

					double t = t_act-period[up_idx]*(time_upstream[up_idx].back()-time_upstream[up_idx][0]);

					if(t >= time_upstream[up_idx][j] && t <= time_upstream[up_idx][j+1])
					{
						got_it=true;
						index_upstream[up_idx] = j;
					}
					else
					{
						j++;
					}
				}

				// interpolating
				double v_h = value_upstream[up_idx][index_upstream[up_idx]+1]; // pressure at higher index
				double v_l = value_upstream[up_idx][index_upstream[up_idx]]; // pressure at lower index
				double t_h = time_upstream[up_idx][index_upstream[up_idx]+1]; // time at higher index
				double t_l = time_upstream[up_idx][index_upstream[up_idx]]; // time at lower index

				double t_in = t_act-period[up_idx]*time_upstream[up_idx].back(); // actual time of the simulation
				double v_in = (v_h-v_l)/(t_h-t_l) * (t_in-t_l) + v_l; // actual pressure of the simulation

				double q_in=0., p_in;
				// there might be several outgoing edge from an upstream node
				for(unsigned int j=0; j<nodes[node_idx[i]]->edge_out.size(); j++)
				{
					int ei = nodes[node_idx[i]]->edge_out[j];
					double dt = edges[ei]->dt_act;
					// calculating vp and pp values of edge
					if(type_upstream[up_idx] == 0) // pressure BC
					{
						p_in = v_in; // v_in is pressure
						q_in += edges[ei]->boundary_pressure_start(dt, p_in);
					}
					else if(type_upstream[up_idx] == 1) // vfr BC
					{
						p_in = edges[ei]->boundary_flowrate_start(dt, v_in);
						q_in += v_in; // v_in is volume flow rate
					}
					else if(type_upstream[up_idx] == 2) // v BC
					{
						double q=0.;
						p_in = edges[ei]->boundary_velocity_start(dt, v_in, q);
						q_in += q;
					}
				}
				for(unsigned int j=0; j<nodes[node_idx[i]]->edge_in.size(); j++)
				{
					int ei = nodes[node_idx[i]]->edge_in[j];
					double dt = edges[ei]->dt_act;
					// calculating vp and pp values of edge
					if(type_upstream[up_idx] == 0) // pressure BC
					{
						p_in = v_in; // v_in is pressure
						q_in += edges[ei]->boundary_pressure_end(dt, p_in);
					}
					else if(type_upstream[up_idx] == 1) // vfr BC
					{
						p_in = edges[ei]->boundary_flowrate_end(dt, v_in);
						q_in += v_in; // v_in is volume flow rate
					}
					else if(type_upstream[up_idx] == 2) // v BC
					{
						double q=0.;
						p_in = edges[ei]->boundary_velocity_end(dt, v_in, q);
						q_in += q;
					}
				}

				if(nodes[node_idx[i]]->do_save_memory)
				{
					nodes[node_idx[i]]->save_field_variables(t_act,p_in,q_in);
				}
			}
			else if(nodes[node_idx[i]]->type_code == 1) // perifera
			{
				double p_out,q_out;
				if(i==0)
				{
					// there can be only one incoming edge
					int edge_index = nodes[node_idx[i]]->edge_out[0];
					double dt = edges[edge_index]->dt_act;
					p_out = nodes[node_idx[i]]->pressure_out;
					q_out = edges[edge_index]->boundary_pressure_start(dt,p_out);
				}
				else
				{
					// there can be only one incoming edge
					int edge_index = nodes[node_idx[i]]->edge_in[0];
					double dt = edges[edge_index]->dt_act;
					p_out = nodes[node_idx[i]]->pressure_out;
					q_out = edges[edge_index]->boundary_pressure_end(dt,p_out);
				}

				if(nodes[node_idx[i]]->do_save_memory)
				{
					nodes[node_idx[i]]->save_field_variables(t_act,p_out,q_out);
				}
			}
			else if(nodes[node_idx[i]]->type_code == 0) // junctions
			{
				// testing if there are concentrated resistances
				bool is_res = false;
				int n1=nodes[node_idx[i]]->edge_in.size();
				int n2=nodes[node_idx[i]]->edge_out.size();

				for(int j=0; j<n1; j++)
				{
					if(edges[nodes[node_idx[i]]->edge_in[j]]->resistance_end != 0.)
					{
						is_res = true;
					}
				}
				for(int j=0; j<n2; j++)
				{
					if(edges[nodes[node_idx[i]]->edge_out[j]]->resistance_start != 0.)
					{
						is_res = true;
					}
				}

				if(!is_res) // NEWEST VERSION with one equation NEWTON (Re=0, Rs=0)
				{
					// initial condition for pressure
					double pp = nodes[node_idx[i]]->pressure.back();

					int k=0;
					double dp=1e10,f;
					do
					{
						// conti equation and derivative
						f = (atmospheric_pressure-pp)/nodes[node_idx[i]]->resistance*nodes[node_idx[i]]->is_resistance;
						double df = -1./nodes[node_idx[i]]->resistance*nodes[node_idx[i]]->is_resistance;

						// edges
						for(int j=0; j<n1; j++)
						{
							int ei = nodes[node_idx[i]]->edge_in[j];
							vector<double> v = edges[ei]->boundary_junction_end(pp, t_act); // [q,dq]
							f  += v[0];
							df += v[1];
						}
						for(int j=0; j<n2; j++)
						{
							int ei = nodes[node_idx[i]]->edge_out[j];
							vector<double> v = edges[ei]->boundary_junction_start(pp, t_act); // [q,dq]
							f  -= v[0];
							df -= v[1];
						}

						dp = -f/df;
						pp += dp;

						k++;
					}
					while(abs(dp) >1e-6 && k<100);

					if(k>=100)
					{
						cout << "\n Newton's technique did NOT converge at node " << nodes[node_idx[i]]->name << " and edge " << edges[e_idx]->name << " for special case (Rs,Re=0)" << endl;
						cout << " k: " << k << " f.norm: " << f << endl;
						//cin.get();
					}

					// substituting back
					nodes[node_idx[i]]->boundary_variables(pp,t_act);

					if(i==1)
					{
						vector<double> v = edges[e_idx]->boundary_junction_end(pp, t_act); // [q,dq]
						double q = v[0];
						edges[e_idx]->boundary_substitute_end(t_act,pp,q);
					}
					else
					{
						vector<double> v = edges[e_idx]->boundary_junction_start(pp, t_act); // [q,dq]
						double q = v[0];
						edges[e_idx]->boundary_substitute_start(t_act,pp,q);
					}
				}
				else// NEW VERSION with NEWTON for GENERAL case
				{
					// initalization
					//edges[e_idx]->y[i](0) = nodes[node_idx[i]]->volume_flow_rate.back()*1.e6; // q [ml/s]
					edges[e_idx]->y[i](0) = nodes[node_idx[i]]->volume_flow_rate.back(); // q [m3/s]
					edges[e_idx]->y[i](1) = nodes[node_idx[i]]->pressure.back(); // p
					for(int j=0; j<n1; j++)
					{
						int ei = nodes[node_idx[i]]->edge_in[j];
						edges[e_idx]->y[i](2+2*j) = edges[ei]->volume_flow_rate_end.back();
						// double A = edges[ei]->area_end.back();
						// edges[e_idx]->y[i](2+2*j+1) = A;
					}

					for(int j=0; j<n2; j++)
					{
						int ei = nodes[node_idx[i]]->edge_out[j];
						edges[e_idx]->y[i](2+2*n1+2*j) = edges[ei]->volume_flow_rate_start.back();
						// double A = edges[ei]->area_start.back();
						// edges[e_idx]->y[i](2+2*n1+2*j+1) = A;
					}

					int n_e_idx;
					int k=0;
					do
					{
						// conti
						edges[e_idx]->f[i](0) = -edges[e_idx]->y[i](0);
						edges[e_idx]->Jac[i](0,0) = -1.;

						// leakage
						edges[e_idx]->f[i](1) = (edges[e_idx]->y[i](1) - atmospheric_pressure)/nodes[node_idx[i]]->resistance - edges[e_idx]->y[i](0); // (p-p0)/R-q
						edges[e_idx]->Jac[i](1,0) = -1.;
						edges[e_idx]->Jac[i](1,1) =  1./nodes[node_idx[i]]->resistance;

						for(int j=0; j<n1; j++) // edges in
						{
							int ei = nodes[node_idx[i]]->edge_in[j];
							// saving edge in index
							if(i==1)
							{
								if(nodes[node_idx[i]]->edge_in[j] == e_idx)
								{
									n_e_idx = j;
								}
							}
							// continouity
							edges[e_idx]->f[i](0) += edges[e_idx]->y[i](2+j);
							edges[e_idx]->Jac[i](0,2+j) = 1.; 

							// char + A
							double qp = edges[e_idx]->y[i](2+j);
							double pp = edges[e_idx]->y[i](1);
							vector<double> v = edges[ei]->boundary_newton_end(qp, pp, t_act); // f_char,dchar_dp,dchard_dq

							// characteristic equation
							edges[e_idx]->f[i](2+2*j) = v[0];
							edges[e_idx]->Jac[i](2+j,1) = v[1]; // dp
							edges[e_idx]->Jac[i](2+j,2+j) = v[2]; // dq
							// edges[e_idx]->Jac[i](2+2*j,2+2*j+1) = -edges[e_idx]->y[i](2+2*j)/(edges[e_idx]->y[i](2+2*j+1)*edges[e_idx]->y[i](2+2*j+1)); // -q/A^2

							// cross-section
							// edges[e_idx]->f[i](2+2*j+1) = v[1];
							// edges[e_idx]->Jac[i](2+2*j+1,1) = v[3]; 
							// edges[e_idx]->Jac[i](2+2*j+1,2+2*j) = v[5]; 
							// edges[e_idx]->Jac[i](2+2*j+1,2+2*j+1) = 1.; 
						}

						for(int j=0; j<n2; j++)
						{
							int ei = nodes[node_idx[i]]->edge_out[j];
							// saving edge out index
							if(i==0)
							{
								if(nodes[node_idx[i]]->edge_out[j] == e_idx)
								{
									n_e_idx = j;
								}
							}

							// continouity
							edges[e_idx]->f[i](0) -= edges[e_idx]->y[i](2+n1+j);
							edges[e_idx]->Jac[i](0,2+2*n1+j) = -1.; 

							// char + A
							double qp = edges[e_idx]->y[i](2+n1+j);
							// double Ap = edges[e_idx]->y[i](2+2*n1+2*j+1);
							double pp = edges[e_idx]->y[i](1);
							vector<double> v = edges[ei]->boundary_newton_start(qp, pp, t_act); // f_char,dchar_dp,dchard_dq

							// characteristic equation
							edges[e_idx]->f[i](2+n1+j) = v[0];
							edges[e_idx]->Jac[i](2+n1+j,1) = v[1]; // dp
							edges[e_idx]->Jac[i](2+n1+j,2+n1+j) = v[2]; // dq
							// edges[e_idx]->Jac[i](2+2*n1+2*j,2+2*n1+2*j+1) = -edges[e_idx]->y[i](2+2*n1+2*j)/(edges[e_idx]->y[i](2+2*n1+2*j+1)*edges[e_idx]->y[i](2+2*n1+2*j+1)); // -q/A^2

							// cross-section
							// edges[e_idx]->f[i](2+2*n1+2*j+1) = v[1];
							// edges[e_idx]->Jac[i](2+2*n1+2*j+1,1) = v[3]; 
							// edges[e_idx]->Jac[i](2+2*n1+2*j+1,2+2*n1+2*j) = v[5];
							// edges[e_idx]->Jac[i](2+2*n1+2*j+1,2+2*n1+2*j+1) = 1.; 
						}

						//cout << "k" << k << endl;
						//cout << "jac " << endl << edges[e_idx]->Jac[i] << endl << endl;
						//cout << "y " << endl << edges[e_idx]->y[i] << endl << endl;
						//cout << "f " << endl << edges[e_idx]->f[i] << endl << endl;
						//cin.get();

						// actually solving Jac*dx = -f
						edges[e_idx]->y[i] += edges[e_idx]->Jac[i].colPivHouseholderQr().solve(-edges[e_idx]->f[i]);

						k++;
					}
					while(edges[e_idx]->f[i].norm() > 1.e-6 && k<100);

					if(k>=100)
					{
						cout << "\n Newton's technique did NOT converge at node " << nodes[node_idx[i]]->name << " and edge " << edges[e_idx]->name << " for general case (Rs,Re!=0)" << endl;
						cout << " k: " << k << " f.norm: " << edges[e_idx]->f[i].norm() << endl;
						//cin.get();
					}

					// saving output vars
					double pp = edges[e_idx]->y[i](1);
					nodes[node_idx[i]]->boundary_variables(pp,t_act);

					if(i==1)
					{
						double q = edges[e_idx]->y[i](2+n_e_idx);
						edges[e_idx]->boundary_substitute_end(t_act,pp,q);
					}
					else
					{
						double q = edges[e_idx]->y[i](2+n1+n_e_idx);
						edges[e_idx]->boundary_substitute_start(t_act,pp,q);
					}
				}
			}
		}
	}
}


//--------------------------------------------------------------
void solver_moc::set_constants(double g, double rho, double nu, double mmHg, double p0, double nu_p, double cfl)
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
		edges[i]->poisson_coefficient = nu_p;
		edges[i]->courant_number = cfl;
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
	for(int j=0; j<type_upstream.size(); j++)
	{
		if(type_upstream[j] == 0)
		{
			for(unsigned int i=0; i<value_upstream[j].size(); i++)
			{
				value_upstream[j][i] *= mmHg_to_Pa;
				value_upstream[j][i] += atmospheric_pressure;
			}
		}
		else if(type_upstream[j] == 1)
		{
			for(unsigned int i=0; i<value_upstream[j].size(); i++)
			{
				value_upstream[j][i] *= 1.e-6;
			}
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



/*//--------------------------------------------------------------
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
}*/

/*//--------------------------------------------------------------
void solver_moc::full_tree()
{
	forward_nodes.clear();
	forward_edges.clear();
	for(int i=0; i<number_of_nodes; i++)
	{
		forward_nodes.push_back(i);
		//if(nodes[i]->is_master_node == false) // not master node
		//{
		//	if(nodes[i]->type_code == 2) // heart
		//	{
		//		nodes[i]->upstream_boundary = true;
		//	}
		//}
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
	//nodes[node_index]->is_upstream_boundary = true;

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
*/

/*//-------------------------------------------------------------- 
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
}*/


/*//--------------------------------------------------------------
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
}*/