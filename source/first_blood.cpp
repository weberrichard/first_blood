#include "first_blood.h"

//--------------------------------------------------------------
first_blood::first_blood(string folder_name)
{
	// saving the name of the folder
	input_folder_path = folder_name;

	// getting rid of global path
	case_name = input_folder_path.substr(input_folder_path.rfind('/')+1);

	// some feedback to console
	cout << " [*] Case name: " << case_name << endl;

	// loading all the input data from csv files
	load_model();

	// setting initial time
	time.push_back(0.);

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
		moc[i]->convert_pt_series();
	}
}

//--------------------------------------------------------------
first_blood::~first_blood(){}

//--------------------------------------------------------------
void first_blood::load_model()
{
	// loading the main csv
	load_main_csv();

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
}

//--------------------------------------------------------------
void first_blood::load_main_csv()
{
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
	}
	else
	{
		cout << "! ERROR !" << endl << " File is not open when calling load_main_csv() function!!! file: " << file_path << "\nExiting..." << endl;
		exit(-1);
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
						moc[j]->load_pt_series(upstream_boundary.file_name);
					}
				}
			}
		}
	}

	file_in.close();
}

//--------------------------------------------------------------
void first_blood::run()
{
	// setting initial conditions
	for(int i=0; i<number_of_moc; i++)
	{
		moc[i]->initialization(pressure_initial);
	}
	for(int i=0; i<number_of_lum; i++)
		
	{
		lum[i]->initialization();
	}

	// setting master nodes
	build_master();

	if(run_type == "forward") // simple forward calculation
	{
		// add everything to forward_*
		for(int i=0; i<number_of_moc; i++)
		{
			moc[i]->full_tree();
		}		

		while(time.back() < time_end)
		{		
			// running mocs
			double dt_master = 1e10;
			for(int i=0; i<number_of_moc; i++)
			{
				double dt = moc[i]->solve_one_step();
				if(dt<dt_master)
				{
					dt_master = dt;
				}
			}

			// if there is no moc we prescribe the dt
			if(number_of_moc == 0)
			{
				dt_master = 1.e-3;
			}

			// basic error handling
			if(dt_master == 1e10)
			{
				cout << " Time step dt is NaN, exiting..." << endl;
				exit(-1);
			}

			// saving time
			time.push_back(time.back()+dt_master);

			// lumped models
			for(int i=0; i<number_of_lum; i++)
			{
				vector<vector<double> > coefs;
				for(int j=0; j<lum[i]->boundary_indices.size(); j++)
				{
					int moc_index = lum[i]->boundary_indices[j][0];
					int moc_edge_index = lum[i]->boundary_indices[j][1];
					int edge_end = lum[i]->boundary_indices[j][2];
					vector<double> c;
					if(edge_end==1)
					{
						c = moc[moc_index]->edges[moc_edge_index]->boundary_master_end(dt_master);
					}
					else
					{
						c = moc[moc_index]->edges[moc_edge_index]->boundary_master_start(dt_master);
					}
					coefs.push_back(c);
				}

				vector<vector<double> > qps = lum[i]->solve_one_step(dt_master, coefs);

				// updating master boundaries in mocs, n. 2
				for(int j=0; j<qps.size(); j++)
				{
					int moc_index = lum[i]->boundary_indices[j][0]; // index of moc
					int moc_edge_index = lum[i]->boundary_indices[j][1]; // edge index in moc
					int edge_end = lum[i]->boundary_indices[j][2];
					int lum_node_index = lum[i]->boundary_indices[j][3]; // node index in lum

					double q = qps[j][0];
					double p = qps[j][1];

					if(edge_end == 1)
					{
						moc[moc_index]->edges[moc_edge_index]->boundary_end_variables(dt_master,p,q);
						int node_index = moc[moc_index]->edges[moc_edge_index]->node_index_end;
						moc[moc_index]->nodes[node_index]->boundary_variables(p);
					}
					else
					{
						moc[moc_index]->edges[moc_edge_index]->boundary_start_variables(dt_master,p,q);
						int node_index = moc[moc_index]->edges[moc_edge_index]->node_index_start;
						moc[moc_index]->nodes[node_index]->boundary_variables(p);
					}
				}
			}

			// postproc: interpolate, save
			for(int i=0; i<number_of_moc; i++)
			{
				moc[i]->post_process(dt_master);
			}			
		}
		
	}
	//else if(run_type == "backward") // backward calculation from inner nodes
	//{
	//}
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
	// LINUX
   mkdir("results",0777);
   mkdir(("results/" + case_name).c_str(),0777);

   // FOR WINDOWS
   //mkdir("results");
   //mkdir(("results/" + folder_name).c_str());

   string folder_name = case_name + "/";

	// saving the results of moc models
	for(int i=0; i<number_of_moc; i++)
	{
		moc[i]->save_results(folder_name + moc[i]->name);
	}

	// saving the results of lumped models
	for(int i=0; i<number_of_lum; i++)
	{
		lum[i]->save_results(folder_name + lum[i]->name);
	}
}


