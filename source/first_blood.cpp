#include "first_blood.h"

//--------------------------------------------------------------
first_blood::first_blood(string folder_name)
{
	// saving the name of the folder
	input_folder_path = folder_name;

	// getting rid of global path
	case_name = input_folder_path.substr(input_folder_path.rfind('/')+1);

	// loading all the input data from csv files
	load_ok = load_model();

	if(load_ok == true)
	{
		// setting constants in every model
		for(int i=0; i<number_of_moc; i++)
		{
			moc[i]->set_constants(gravity, density, kinematic_viscosity, mmHg_to_Pa, atmospheric_pressure,poisson_coefficient, courant_number);
		}

		for(int i=0; i<number_of_lum; i++)
		{
			lum[i]->set_constants(gravity, density, kinematic_viscosity, mmHg_to_Pa, atmospheric_pressure);
		}
		// converting moc t-p to SI
		for(int i=0; i<number_of_moc; i++)
		{
			moc[i]->convert_time_series();
		}
	}
}

//--------------------------------------------------------------
first_blood::~first_blood()
{
	initialization();
}

//--------------------------------------------------------------
bool first_blood::load_model()
{
	// loading the main csv
	bool load_ok = load_main_csv();
	if(load_ok == false)
	{
		return load_ok;
	}

	// setting the number of models
	number_of_moc = moc.size();
	number_of_lum = lum.size();
	number_of_nodes = nodes.size();

	// loading the csv for moc models
	for(int i=0; i<number_of_moc; i++)
	{
		moc[i]->load_model();
	}

	// loading the csv for lumped models
	for(int i=0; i<number_of_lum; i++)
	{
		lum[i]->load_model();
	}
	return load_ok;
}

//--------------------------------------------------------------
bool first_blood::load_main_csv()
{
	load_ok = false;
	ifstream file_in;
	string file_path = input_folder_path + "/main.csv";
	file_in.open(file_path);
	string line;
	if(file_in.is_open())
	{
		int nm=0,nl=0,nn=0; // nm for moc, nl for lumped, nn for nodes
		while(getline(file_in,line))
		{
			// clearing spaces
			line.erase(remove(line.begin(), line.end(), ' '), line.end());
			line.erase(remove(line.begin(), line.end(), '\n'), line.end());
			line.erase(remove(line.begin(), line.end(), '\r'), line.end());

			// separating strings to vector by comma
			vector<string> sv = separate_line(line);

			if(sv[0] == "run") // run type, boundaries if given
			{
				run_type = sv[1];
				if(sv.size()>2)
				{
					upstream_boundary.node = sv[2];
					upstream_boundary.file_name = sv[3];
				}
				else
				{
					upstream_boundary.node = "";
					upstream_boundary.file_name = "";
				}
			}
			else if(sv[0] == "time") // time
			{
				time_end = stod(sv[1],0);
			}
			else if(sv[0] == "material") // material model
			{
				if(sv[1] == "linear")
				{
					material_type = 0;
				}
				else if(sv[1] == "Olufsen" || sv[1] == "olufsen")
				{
					material_type = 1;
				}
			}
			else if(sv[0] == "solver") // setting solver type, 0: maccormack, 1: moc
			{
				if(sv[1] == "maccormack" || sv[1] == "Maccormack" || sv[1] == "MacCormack")
				{
					solver_type = 0;
				}
				else if(sv[1] == "moc")
				{
					solver_type = 1;
				}
				else
				{
					cout << "Solver type: " << sv[1] << " is not known. Using MacCormack, continoiug..." << endl;
				}
			}
			else if(sv[0] == "moc") // 1D moc model
			{
				moc.push_back(new solver_moc(sv[1],input_folder_path));
				int k=2;
				while(k<sv.size())
				{
					moc[nm]->boundary_main_node.push_back(sv[k]);
					moc[nm]->boundary_model_node.push_back(sv[k+1]);
					k+=2;
				}
				nm++;
			}
			else if(sv[0] == "lumped" || sv[0] == "lum") // 0D lumed models
			{
				lum.push_back(new solver_lumped(sv[1],input_folder_path));
				int k=2;
				while(k<sv.size())
				{
					lum[nl]->boundary_main_node.push_back(sv[k]);
					lum[nl]->boundary_model_node.push_back(sv[k+1]);
					k+=2;
				}
				nl++;
			}
			else if(sv[0] == "node") // main nodes between models
			{
				nodes.push_back(sv[1]);
				nn++;
			}
		}
		load_ok = true;
	}
	else
	{
		cout << "! ERROR !" << endl << " File is not open when calling load_main_csv() function!!! file: " << file_path << "\n" << endl;
		load_ok = false;
		return load_ok;
	}

	// loading p-t if it is given in main.csv
	if(upstream_boundary.node != "")
	{
		for(int i=0; i<nodes.size(); i++)
		{
			for(int j=0; j<moc.size(); j++)
			{
				for(int k=0; k<moc[j]->nodes.size(); k++)
				{
					if(nodes[i] == moc[j]->nodes[k]->name)
					{
						moc[j]->load_time_series(upstream_boundary.file_name);
					}
				}
			}
		}
	}

	file_in.close();

	return load_ok;
}

//--------------------------------------------------------------
bool first_blood::run()
{
	bool is_run_ok;

	// initialization of the whole model with mocs and lums
	initialization();

	if(init_from_file)
	{
		load_initials();
	}

	if(run_type == "forward") // simple forward calculation
	{
		is_run_ok = true;

		if(number_of_moc>0)
		{
			// initial timesteps
			int moc_idx=0, e_idx;
			for(int i=0; i<number_of_moc; i++)
			{
				moc[i]->timesteps();
			}	

			// finding lowest new timestep
			double t_act = lowest_new_time(moc_idx, e_idx);
			double t_old = -1.e10;


			// main cycle
			while(!is_run_end(t_act,t_old) && is_run_ok)
			{
				
				//cout << "t: " << t_act << endl;
				/*if(t_act>7.) // improve this
				{
					int idx = lum_id_to_index("heart_kim");
					heart_rate = 90.;
					lum[idx]->heart_rate = heart_rate;
					time_period = 60./heart_rate;
				}

				if(t_act>13.)
				{
					do_autoregulation = true;
				}

				if(do_autoregulation)
				{
					autoregulation();
				}*/

				// solving lowest edge inner points
				if(solver_type == 0)
				{
					moc[moc_idx]->edges[e_idx]->solve_maccormack();
				}
				else if(solver_type == 1)
				{
					moc[moc_idx]->edges[e_idx]->solve_moc();
				}

				// boundaries (inner)
				moc[moc_idx]->boundaries(e_idx, t_act);

				// boundaries with 0D if exist
				int si = moc[moc_idx]->edges[e_idx]->node_index_start;
				if(moc[moc_idx]->nodes[si]->is_master_node)
				{
					int lum_idx = moc[moc_idx]->nodes[si]->master_node_lum;
					solve_lum_newton(lum_idx, t_act);
				}

				int ei = moc[moc_idx]->edges[e_idx]->node_index_end;
				if(moc[moc_idx]->nodes[ei]->is_master_node)
				{
					int lum_idx = moc[moc_idx]->nodes[ei]->master_node_lum;
					solve_lum_newton(lum_idx, t_act);
				}

				// postproc: interpolate, save
				moc[moc_idx]->edges[e_idx]->update();
				//moc[moc_idx]->edges[e_idx]->update_variables();
				if(moc[moc_idx]->edges[e_idx]->do_save_memory)
				{
					moc[moc_idx]->edges[e_idx]->save_field_variables();
				}

				// get the time average values, e.g MAP
				//if(moc[moc_idx]->edges[e_idx]->ID == "A1")
				//{
				//	calculate_time_average();
				//}

				// new timestep
				is_run_ok = moc[moc_idx]->edges[e_idx]->new_timestep();

				// find new lowest timestep
				t_old = t_act;
				t_act = lowest_new_time(moc_idx, e_idx);

				// decreasing time_counter if it passed t_end
				if(t_act >= time_end)
				{
					time_counter--;
				}

				// update the period
				period = floor(t_act/time_period);
			}
		}
		else // only LUMPED MODEL without any moc
		{
			if(number_of_lum==0)
			{
				cout << " There is no LUMPED model (nor MOC).\n Exiting..." << endl;
				return false;
			}

			double t_act = 0.;
			while(t_act<time_end)
			{
				t_act = lum[0]->time.back() + dt_lumped;
				solve_lum_newton(0, t_act);
			}
		}
	}

	return is_run_ok;
}

//--------------------------------------------------------------
double first_blood::lowest_new_time(int &moc_idx, int &e_idx)
{
	double t_act=1.e10;
	int idx;
	for(int i=0; i<number_of_moc; i++)
	{
		double t = moc[i]->min_time(idx);
		if(t<t_act)
		{
			t_act = t;
			moc_idx = i;
			e_idx = idx;
		}
	}

	return t_act;
}

//--------------------------------------------------------------
void first_blood::solve_lum_newton(int index, double t_act)
{
	/*
	x = [q1,q2,...qm,p1,p2,p3,...pn,y1,y2,...ye,qmoc1,...qmock] y: for elastance if present
	f = [edge1,edge2,...edgem,node1,node2,...node,elas1,elas2,...elase,char1,char2,...]
	*/
	int m = lum[index]->number_of_edges;
	int n = lum[index]->number_of_nodes;
	int l = lum[index]->number_of_elastance;
	int k = lum[index]->boundary_indices.size();

	// setting initial conditions for lumped part
	lum[index]->initialization_newton(t_act);
	// setting initial conditions for moc part
	for(int j=0; j<lum[index]->boundary_indices.size(); j++)
	{
		int moc_index = lum[index]->boundary_indices[j][0];
		int moc_edge_index = lum[index]->boundary_indices[j][1];
		int edge_end = lum[index]->boundary_indices[j][2];
		int N = m+n+2*l+j;
		moc[moc_index]->initialization_newton(lum[index]->x,N,moc_edge_index,edge_end);
	}

	// start of newton iteration
	int i=0;
	do
	{
		lum[index]->coefficients_newton(t_act);

		for(int j=0; j<lum[index]->boundary_indices.size(); j++)
		{
			int moc_index = lum[index]->boundary_indices[j][0];
			int moc_edge_index = lum[index]->boundary_indices[j][1];
			int edge_end = lum[index]->boundary_indices[j][2];
			int lum_node_index = lum[index]->boundary_indices[j][3];
			double q = lum[index]->x(n+m+2*l+j);
			// double A = lum[index]->x(n+m+2*l+2*j+1);
			double p = lum[index]->x(m+lum_node_index)*mmHg_to_Pa;
			vector<double> v; // f_char,dchar_dp,dchard_dq
			if(edge_end==1)
			{
				v = moc[moc_index]->edges[moc_edge_index]->boundary_newton_end(q,p,t_act);
			}
			else
			{
				v = moc[moc_index]->edges[moc_edge_index]->boundary_newton_start(q,p,t_act);
			}

			// node continouity equation and jacobian
			lum[index]->f(m+lum_node_index) += edge_end*lum[index]->x(m+n+2*l+j)*1.e6; // sign(q)*q*1e6 [ml/s]
			lum[index]->Jac(m+lum_node_index,m+n+2*l+j) = edge_end*1.e6; // for node continouity // sign(q)*1e6 [1]

			// characteristic equation
			// also converting to non_SI for numerical stability
			lum[index]->f(m+n+2*l+j) = v[0];
			lum[index]->Jac(m+n+2*l+j,m+lum_node_index) = v[1]*mmHg_to_Pa; // dp
			lum[index]->Jac(m+n+2*l+j,m+n+2*l+j) = v[2]; // dQ
			// lum[index]->Jac(m+n+2*l+2*j,m+n+2*l+2*j+1) = -q/A/A; // dA

			// cross section - pressure equation
			// lum[index]->f(m+n+2*l+2*j+1) = v[1];
			// lum[index]->Jac(m+n+2*l+2*j+1,m+lum_node_index) = v[3]*mmHg_to_Pa; // dp
			// lum[index]->Jac(m+n+2*l+2*j+1,m+n+2*l+2*j) = v[5]; // dQ
			// lum[index]->Jac(m+n+2*l+2*j+1,m+n+2*l+2*j+1) = 1.; // dA
		}

		// cout << endl << i << " ITERATION " << i << endl;
		// cout << endl << lum[index]->Jac << endl;
		// cout << "x" << endl << lum[index]->x << endl;
		// cout << "f" << endl << lum[index]->f << endl;

		// actually solving the Newton's technique
		VectorXd dx = lum[index]->Jac.colPivHouseholderQr().solve(-lum[index]->f);
		lum[index]->x += dx;

		//cout << " ITERATION END " << endl;
		//cout << endl;
		//cin.get();
		
		i++;
	}
	while(lum[index]->f.norm() > 1e-5 && i<100);

	if(i>=100)
	{
		cout << "\n !!! ERROR !!! Newton's technique did NOT converge at Lum model " << lum[index]->name << endl;
		cout << endl << i << " ITERATION " << i << endl;
		cout << lum[index]->Jac << endl;
		cout << "x" << endl << lum[index]->x << endl;
		cout << "f" << endl << lum[index]->f << endl;
		exit(-1);
	}

	// substitute the results back to 0D
	lum[index]->substitute_newton(t_act);

	// substitute the results back to 1D
	for(int j=0; j<lum[index]->boundary_indices.size(); j++)
	{
		int moc_index = lum[index]->boundary_indices[j][0];
		int moc_edge_index = lum[index]->boundary_indices[j][1];
		int edge_end = lum[index]->boundary_indices[j][2];
		int lum_node_index = lum[index]->boundary_indices[j][3];
		double p = lum[index]->x(m+lum_node_index)*mmHg_to_Pa;
		double q = lum[index]->x(m+n+2*l+j);
		moc[moc_index]->substitute_newton(moc_edge_index,edge_end,t_act,p,q);
	}
}

//--------------------------------------------------------------
void first_blood::calculate_time_average()
{
	// mean arterial pressure
	string id = "A1";
	int idx = moc[0]->edge_id_to_index(id);
	double tn = moc[0]->edges[idx]->time.back();
	double vn = moc[0]->edges[idx]->pressure_start.back();
	map->update(tn,vn,time_period);

	// mean cerebral flow rate
	vector<string> ids{"A5","A6","A20","A15"};
	vn=0.;
	for(int i=0; i<ids.size(); i++)
	{
		idx = moc[0]->edge_id_to_index(ids[i]);
		vn += moc[0]->edges[idx]->volume_flow_rate_start.back();
	}
	cfr->update(tn,vn,time_period);
}

//--------------------------------------------------------------
void first_blood::autoregulation()
{	
	//double S = 1.;
	//double pset = 90.;
	//pset = pset*mmHg_to_Pa + atmospheric_pressure;
	//double factor = 1.0 + S*(map->average.back()-pset);
	//cout << "f: " << factor << endl;

	double pset = 90.;
	double factor = (map->average.back()-atmospheric_pressure)/mmHg_to_Pa / pset;

	vector<string> perif_brain{"p25","p26","p27","p28","p29","p30","p31","p32","p33","p34","p35","p36","p37","p38","p39","p40","p41","p42","p43","p44","p45","p46"};

	for(int i=0; i<perif_brain.size(); i++)
      {
         string id = perif_brain[i];
         int idx = lum_id_to_index(id);
         lum[idx]->edges[0]->parameter_factor = factor;
         lum[idx]->edges[1]->parameter_factor = factor;
         lum[idx]->edges[2]->parameter_factor = 1./factor;
      }
}

//--------------------------------------------------------------
void first_blood::save_time_average(string folder_name)
{
	mkdir("results",0777);
	mkdir(("results/"+folder_name).c_str(),0777);

	string file_name = folder_name + "/arterial/map.txt";
	map->save_results(file_name);

	file_name = folder_name + "/arterial/cfr.txt";
	cfr->save_results(file_name);
}

//--------------------------------------------------------------
void first_blood::save_time_average(double dt, string folder_name)
{
	mkdir("results",0777);
	mkdir(("results/"+folder_name).c_str(),0777);

	string file_name = folder_name + "/arterial/map.txt";
	map->save_results(dt, file_name);

	file_name = folder_name + "/arterial/cfr.txt";
	cfr->save_results(dt, file_name);
}

//--------------------------------------------------------------
bool first_blood::is_run_end(double t_act, double t_old)
{
	if(is_periodic_run)
	{
		if(t_act<time_end_min)
		{
			return false;
		}
		else if(t_act>time_end_max)
		{
			time_end = time_end_max;
			return true;
		}
		else
		{
			int n=0;
			double t = t_act;
			while(t >= time_period)
			{
				t -= time_period;
				n++;
			}

 			// checking the end of a cycle
 			if((double)n*time_period>=t_old && (double)n*time_period<t_act)
 			{
				int idx = moc[0]->node_id_to_index(time_node);
				if(idx<0)
				{
					time_node = moc[0]->nodes[0]->name;
					idx = 0;
				}
				if(time_var == "P")
				{
					double val = (systole(moc[0]->nodes[idx]->pressure, moc[0]->nodes[idx]->time, (n-1.)*time_period)-atmospheric_pressure)/mmHg_to_Pa;
					if(abs((val-time_val_old)/time_val_old)<.0001)
					{
						time_end = t_act;
						return true;
					}
					else
					{
						time_val_old = val;
						return false;
					}
				}
				else
				{
					cout << "\n time_var: " << time_var << " is not valid, avaialble: P" << endl;
					time_end = t_act;
					return true;
				}
 			}
 			else
 			{
 				return false;
 			}
		}
	}
	else
	{
		if(time_counter>0)
		{
			return false;
		}
		else
		{
			return true;
		}
	}
}

//--------------------------------------------------------------
void first_blood::initialization()
{
	// fb class
	number_of_moc = moc.size();
	number_of_lum = lum.size();
	number_of_nodes = nodes.size();
	time_counter = 0;

	// setting initial conditions
	for(int i=0; i<number_of_moc; i++)
	{
		moc[i]->initialization(pressure_initial,material_type);
		time_counter += moc[i]->number_of_edges;
	}
	for(int i=0; i<number_of_lum; i++)
	{
		lum[i]->initialization(heart_rate);
		time_counter += 1;
	}

	// setting master nodes
	build_master();

	// setting number of moc in lum
	for(int i=0; i<number_of_lum; i++)
	{
		lum[i]->set_newton_size();
	}

	// making sure time_node is saved
	if(is_periodic_run)
	{
		int idx = moc[0]->node_id_to_index(time_node);
		if(idx<0)
		{
			time_node = moc[0]->nodes[0]->name;
			idx = 0;
		}
		if(moc.size()>0)
		{
			vector<string> el, nl{time_node};
		   set_save_memory(moc[0]->name,"moc",el,nl);
		}
	}

	// saving these for time average vectors
	// vector<string> el{"A1","A5","A6","A15","A20"}, nl;
	// set_save_memory(moc[0]->name,"moc",el,nl);

	// time average stuff
	map = new time_average();
	cfr = new time_average();
}

//--------------------------------------------------------------
void first_blood::build_master()
{
	for(int i=0; i<nodes.size(); i++)
	{
		for(int j=0; j<number_of_lum; j++)
		{
			for(int k=0; k<lum[j]->boundary_main_node.size(); k++)
			{
				if(nodes[i] == lum[j]->boundary_main_node[k])
				{
					int kk = lum[j]->node_id_to_index(lum[j]->boundary_model_node[k]);
					lum[j]->nodes[kk]->is_master_node = true;

					for(int l=0; l<number_of_moc; l++)
					{
						for(int m=0; m<moc[l]->boundary_main_node.size(); m++)
						{
							if(nodes[i] == moc[l]->boundary_main_node[m])
							{
								int mm = moc[l]->node_id_to_index(moc[l]->boundary_model_node[m]);
								moc[l]->nodes[mm]->is_master_node = true;
								moc[l]->nodes[mm]->master_node_lum = j;

								for(int n=0; n<moc[l]->nodes[mm]->edge_in.size(); n++)
								{
									vector<int> v{l,moc[l]->nodes[mm]->edge_in[n],1,kk};
									// index of moc, index of edge in moc, start(-1)/end(+1), index of node in lumped
									lum[j]->boundary_indices.push_back(v);
								}
								for(int n=0; n<moc[l]->nodes[mm]->edge_out.size(); n++)
								{
									vector<int> v{l,moc[l]->nodes[mm]->edge_out[n],-1,kk};
									// index of moc, index of edge in moc, start(-1)/end(+1), index of node in lumped
									lum[j]->boundary_indices.push_back(v);
								}
							}
						}
					}
				}
			}
		}
	}
}

//--------------------------------------------------------------
void first_blood::save_results()
{
   mkdir("results",0777);
   mkdir(("results/" + case_name).c_str(),0777);

   string folder_name = case_name;

	// saving the results of moc models
	for(int i=0; i<number_of_moc; i++)
	{
		moc[i]->save_results(folder_name);
	}

	// saving the results of lumped models
	for(int i=0; i<number_of_lum; i++)
	{
		lum[i]->save_results(folder_name);
	}

	// saving time averages, e.g. map, cfr
	// save_time_average("results/" + folder_name);
}

//--------------------------------------------------------------
void first_blood::save_results(string folder_name)
{
   mkdir("results",0777);
   mkdir(("results/" + folder_name).c_str(),0777);

	// saving the results of moc models
	for(int i=0; i<number_of_moc; i++)
	{
		moc[i]->save_results(folder_name);
	}

	// saving the results of lumped models
	for(int i=0; i<number_of_lum; i++)
	{
		lum[i]->save_results(folder_name);
	}

	// saving time averages, e.g. map, cfr
	// save_time_average("results/" + folder_name);
}

//--------------------------------------------------------------
void first_blood::save_results(double dt)
{
   mkdir("results",0777);
   mkdir(("results/" + case_name).c_str(),0777);

   string folder_name = case_name;

	// saving the results of moc models
	for(int i=0; i<number_of_moc; i++)
	{
		moc[i]->save_results(dt, folder_name);
	}

	// saving the results of lumped models
	for(int i=0; i<number_of_lum; i++)
	{
		lum[i]->save_results(dt, folder_name);
	}

	// saving time averages, e.g. map, cfr
	save_time_average(dt, "results/" + folder_name);
}

//--------------------------------------------------------------
void first_blood::save_results(double dt, string folder_name)
{
   mkdir("results",0777);
   mkdir(("results/" + folder_name).c_str(),0777);

	// saving the results of moc models
	for(int i=0; i<number_of_moc; i++)
	{
		moc[i]->save_results(dt, folder_name);
	}

	// saving the results of lumped models
	for(int i=0; i<number_of_lum; i++)
	{
		lum[i]->save_results(dt, folder_name);
	}
}

//--------------------------------------------------------------
void first_blood::save_results(string folder_name, string model_name, string model_type, vector<string> edge_list, vector<string> node_list)
{
   mkdir("results",0777);
   mkdir(("results/" + folder_name).c_str(),0777);

	if(model_type == "moc")
	{
		for(int i=0; i<moc.size(); i++)
		{
			if(model_name == moc[i]->name)
			{
				moc[i]->save_results(folder_name,edge_list,node_list);
			}
		}
	}
	else if(model_type == "lum")
	{
		for(int i=0; i<lum.size(); i++)
		{
			if(model_name == lum[i]->name)
			{
				lum[i]->save_results(folder_name,edge_list,node_list);			
			}
		}
	}
}

//--------------------------------------------------------------
void first_blood::save_results(double dt, string folder_name, string model_name, string model_type, vector<string> edge_list, vector<string> node_list)
{
	mkdir("results",0777);
   mkdir(("results/" + folder_name).c_str(),0777);

	if(model_type == "moc")
	{
		for(int i=0; i<moc.size(); i++)
		{
			if(model_name == moc[i]->name)
			{
				moc[i]->save_results(dt,folder_name,edge_list,node_list);
			}
		}
	}
	else if(model_type == "lum")
	{
		for(int i=0; i<lum.size(); i++)
		{
			if(model_name == lum[i]->name)
			{
				lum[i]->save_results(dt,folder_name,edge_list,node_list);			
			}
		}
	}
}

//--------------------------------------------------------------
void first_blood::save_model(string model_name)
{
	string folder_name = "../../models/";

	save_model(model_name, folder_name);
}

//--------------------------------------------------------------
void first_blood::save_model(string model_name, string folder_name)
{
   mkdir(folder_name.c_str(),0777);
   mkdir((folder_name+"/"+model_name).c_str(),0777);

	// saving moc models
	for(int i=0; i<number_of_moc; i++)
	{
		moc[i]->save_model(model_name, folder_name);
	}

	// saving lumped models
	for(int i=0; i<number_of_lum; i++)
	{
		lum[i]->save_model(model_name, folder_name);
	}

	// saving main model TODO
	FILE *out_file;
	string file_name = folder_name + "/" + model_name + "/main.csv";
	out_file = fopen(file_name.c_str(),"w");

	fprintf(out_file, "run,forward\n");
	fprintf(out_file, "time,%6.3f\n",time_end);
	fprintf(out_file, "\n");

	fprintf(out_file, "type,name,main node,model node,main node,model node,...\n");
	for(int i=0; i<number_of_moc; i++)
	{
		fprintf(out_file, "moc,%s",moc[i]->name.c_str());
		for(int j=0; j<moc[i]->boundary_main_node.size(); j++)
		{
			fprintf(out_file, ",%s,%s",moc[i]->boundary_main_node[j].c_str(),moc[i]->boundary_model_node[j].c_str());
		}
		fprintf(out_file, "\n");
	}
	fprintf(out_file, "\n");

	for(int i=0; i<number_of_lum; i++)
	{
		fprintf(out_file, "lumped,%s",lum[i]->name.c_str());
		for(int j=0; j<lum[i]->boundary_main_node.size(); j++)
		{
			fprintf(out_file, ",%s,%s",lum[i]->boundary_main_node[j].c_str(),lum[i]->boundary_model_node[j].c_str());
		}
		fprintf(out_file, "\n");
	}
	fprintf(out_file, "\n");

	for(int i=0; i<number_of_nodes; i++)
	{
		fprintf(out_file, "node,%s\n",nodes[i].c_str());
	}
	fprintf(out_file, "\n");
   fclose(out_file);
}

//--------------------------------------------------------------
void first_blood::clear_save_memory()
{
	for(int i=0; i<number_of_moc; i++)
	{
		moc[i]->clear_save_memory();
	}
	for(int i=0; i<number_of_lum; i++)
	{
		lum[i]->clear_save_memory();
	}
}

//--------------------------------------------------------------
void first_blood::set_save_memory(string model_name, string model_type, vector<string> edge_list, vector<string> node_list)
{
	if(model_type == "moc")
	{
		for(int i=0; i<number_of_moc; i++)
		{
			if(moc[i]->name == model_name)
			{
				moc[i]->set_save_memory(edge_list, node_list);
			}
		}
	}
	if(model_type == "lum" || model_type == "lumped")
	{
		for(int i=0; i<number_of_lum; i++)
		{
			if(lum[i]->name == model_name)
			{
				lum[i]->set_save_memory(edge_list, node_list);
			}
		}
	}
}

//--------------------------------------------------------------
int first_blood::lum_id_to_index(string lum_id)
{
	int i=0, idx=-1;
	bool got_it=false;
	while(i<number_of_lum && !got_it)
	{
		if(lum_id.compare(lum[i]->name) == 0)
		{
			got_it = true;
			idx = i;
		}
		i++;
	}
	if(idx == -1)
	{
		cout << "\n !!!WARNING!!!\n solver_moc::lum_id_to_index function\nLum model is not existing, lum_id: " << lum_id << "\n Continouing..." << endl;
	}
	return idx;
}

//--------------------------------------------------------------
void first_blood::save_initials(string model_name, string folder_name)
{
   mkdir(folder_name.c_str(),0777);
   mkdir((folder_name+"/"+model_name).c_str(),0777);
   mkdir((folder_name+"/"+model_name+"/init").c_str(),0777);

   for(int i=0; i<number_of_moc; i++)
   {
   	moc[i]->save_initials(model_name, folder_name);
   }

   for(int i=0; i<number_of_lum; i++)
   {
   	lum[i]->save_initials(model_name, folder_name);
   }
}

//--------------------------------------------------------------
void first_blood::load_initials()
{
	for(int i=0; i<number_of_moc; i++)
	{
		moc[i]->load_initials();
	}

	for(int i=0; i<number_of_lum; i++)
	{
		lum[i]->load_initials();
	}
}

//--------------------------------------------------------------
/*void first_blood::solve_lum(int index, double dt)
{
	vector<vector<double> > coefs;
	for(int j=0; j<lum[index]->boundary_indices.size(); j++)
	{
		int moc_index = lum[index]->boundary_indices[j][0];
		int moc_edge_index = lum[index]->boundary_indices[j][1];
		int edge_end = lum[index]->boundary_indices[j][2];
		vector<double> c;
		if(edge_end==1)
		{
			c = moc[moc_index]->edges[moc_edge_index]->boundary_master_end(dt);
		}
		else
		{
			c = moc[moc_index]->edges[moc_edge_index]->boundary_master_start(dt);
		}
		coefs.push_back(c);
	}

	vector<vector<double> > qps = lum[index]->solve_one_step(dt, coefs);

	// updating master boundaries in mocs, n. 2
	for(int j=0; j<qps.size(); j++)
	{
		int moc_index = lum[index]->boundary_indices[j][0]; // index of moc
		int moc_edge_index = lum[index]->boundary_indices[j][1]; // edge index in moc
		int edge_end = lum[index]->boundary_indices[j][2];
		int lum_node_index = lum[index]->boundary_indices[j][3]; // node index in lum

		double q = qps[j][0];
		double p = qps[j][1];

		double t_act = moc[moc_index]->edges[moc_edge_index]->time.back() + dt;
		if(edge_end == 1)
		{
			moc[moc_index]->edges[moc_edge_index]->boundary_end_variables(dt,p,q);
			int node_index = moc[moc_index]->edges[moc_edge_index]->node_index_end;
			moc[moc_index]->nodes[node_index]->boundary_variables(p,t_act);
		}
		else
		{
			moc[moc_index]->edges[moc_edge_index]->boundary_start_variables(dt,p,q);
			int node_index = moc[moc_index]->edges[moc_edge_index]->node_index_start;
			moc[moc_index]->nodes[node_index]->boundary_variables(p,t_act);
		}
	}
}*/
