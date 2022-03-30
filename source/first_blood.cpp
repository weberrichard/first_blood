#include "first_blood.h"

//--------------------------------------------------------------
first_blood::first_blood(string folder_name)
{
	// saving the name of the folder
	input_folder_path = folder_name;

	// getting rid of global path
	case_name = input_folder_path.substr(input_folder_path.rfind('/')+1);

	// some feedback to console
	//cout << " [*] Case name: " << case_name << endl;

	// loading all the input data from csv files
	load_ok = load_model();

	if(load_ok == true)
	{
		// setting constants in every model
		for(int i=0; i<number_of_moc; i++)
		{
			moc[i]->set_constants(gravity, density, kinematic_viscosity, mmHg_to_Pa, atmospheric_pressure, beta);
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
			else if(sv[0] == "lumped") // 0D lumed models
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
		//cout << "! ERROR !" << endl << " File is not open when calling load_main_csv() function!!! file: " << file_path << "\n" << endl;
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

	if(run_type == "forward") // simple forward calculation
	{
		is_run_ok = true;

		// add everything to forward_*
		for(int i=0; i<number_of_moc; i++)
		{
			moc[i]->full_tree();
		}

		// initial timesteps
		int moc_idx=0, e_idx;
		for(int i=0; i<number_of_moc; i++)
		{
			double t = moc[i]->timesteps(e_idx);
		}
		double t_act = lowest_new_time(moc_idx, e_idx);


		// main cycle
		while(time_counter > 0)
		{
			// solving lowest edge inner points
			moc[moc_idx]->edges[e_idx]->solve();

			// boundaries (inner + lumped if it is)
			moc[moc_idx]->boundaries(e_idx, t_act);

			int si = moc[moc_idx]->edges[e_idx]->node_index_start;
			if(moc[moc_idx]->nodes[si]->is_master_node)
			{
				double dt = moc[moc_idx]->edges[e_idx]->dt_act;
				int lum_idx = moc[moc_idx]->nodes[si]->master_node_lum;
				solve_lum(lum_idx, dt);
			}

			int ei = moc[moc_idx]->edges[e_idx]->node_index_end;
			if(moc[moc_idx]->nodes[ei]->is_master_node)
			{
				double dt = moc[moc_idx]->edges[e_idx]->dt_act;
				int lum_idx = moc[moc_idx]->nodes[ei]->master_node_lum;
				solve_lum(lum_idx, dt);
			}

			// postproc: interpolate, save
			moc[moc_idx]->edges[e_idx]->interpolate();
			moc[moc_idx]->edges[e_idx]->update_variables();
			if(moc[moc_idx]->edges[e_idx]->do_save_memory)
			{
				moc[moc_idx]->edges[e_idx]->save_field_variables();
			}

			// new timestep
			moc[moc_idx]->edges[e_idx]->new_timestep();

			// find new lowest timestep
			t_act = lowest_new_time(moc_idx, e_idx);

			// decreasing time_counter if it passed t_end
			if(t_act >= time_end)
			{
				time_counter--;
			}
		}
	}

	return is_run_ok;
}

//--------------------------------------------------------------
void first_blood::solve_lum(int index, double dt)
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
}

//--------------------------------------------------------------
double first_blood::lowest_new_time(int &moc_idx, int &e_idx)
{
	double tact=1.e10;
	for(int i=0; i<number_of_moc; i++)
	{
		double t = moc[i]->min_time(e_idx);
		if(t<tact)
		{
			tact = t;
			moc_idx = i;
		}
	}

	return tact;
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
		moc[i]->initialization(pressure_initial);
		time_counter += moc[i]->number_of_edges;
	}
	for(int i=0; i<number_of_lum; i++)
	{
		lum[i]->initialization();
		time_counter += 1;
	}

	// setting master nodes
	build_master();
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

								int n;
								int nn;
								if(moc[l]->nodes[mm]->edge_in.size() != 0)
								{
									n = moc[l]->nodes[mm]->edge_in[0];
									nn = 1;
								}
								else
								{
									n = moc[l]->nodes[mm]->edge_out[0];
									nn = -1;
								}
								// index of moc, index of edge in moc, start(-1)/end(+1), index of node in lumped
								vector<int> v{l,n,nn,kk};
								lum[j]->boundary_indices.push_back(v);
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
}

//--------------------------------------------------------------
void first_blood::save_results(string folder_name, string model_name, string model_type, vector<string> edge_list, vector<string> node_list)
{
   mkdir("results",0777);
   mkdir(("results/" + case_name).c_str(),0777);

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
   mkdir(("results/" + case_name).c_str(),0777);

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
   mkdir((folder_name+model_name).c_str(),0777);

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
	string file_name = folder_name + model_name + "/main.csv";
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
