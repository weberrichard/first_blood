#include "solver_moc.h"

//--------------------------------------------------------------
void solver_moc::load_model()
{
	ifstream file_in;
	string file_name = input_folder_path + '/' + name + ".csv";
	file_in.open(file_name);
	string line;
	if(file_in.is_open())
	{
		int i=0,j=0; // i for edges, j for nodes
		while(getline(file_in,line))
		{
			// clearing spaces and \n
			line.erase(remove(line.begin(), line.end(), ' '), line.end());
			line.erase(remove(line.begin(), line.end(), '\n'), line.end());
			line.erase(remove(line.begin(), line.end(), '\r'), line.end());

			// seperating the strings by comma
			vector<string> sv = separate_line(line);
			
			if(sv[0] == "vis" || sv[0] == "visM") // edges
			{
				edges.push_back(new moc_edge(sv[1]));
				edges[i]->name = sv[2];
				edges[i]->node_name_start = sv[3];
				edges[i]->node_name_end = sv[4];
				edges[i]->nominal_diameter_start = stod(sv[5],0);
				edges[i]->nominal_diameter_end = stod(sv[6],0);
				edges[i]->nominal_thickness_start = stod(sv[7],0);
				edges[i]->nominal_thickness_end = stod(sv[8],0);
				edges[i]->length = stod(sv[9],0);
				edges[i]->division_points = stoi(sv[10],0);
				edges[i]->elasticity_spring = stod(sv[11],0);
				edges[i]->elasticity_voigt = stod(sv[12],0);
				edges[i]->viscosity = stod(sv[13],0);
				edges[i]->resistance_start = stod(sv[14],0);
				edges[i]->resistance_end = stod(sv[15],0);
				i++;
			}
			else if(sv[0] == "elag" || sv[0] == "junction" || sv[0] == "node") // nodes
			{
				nodes.push_back(new moc_node(sv[1]));
				if(sv[3] == "")
				{
					nodes[j]->resistance = 1.; // this will be multiplied with is_resistance i.e. there is no leakage
					nodes[j]->is_resistance = 0.;
				}
				else
				{
					double R = stod(sv[3],0);
					if(R == 0)
					{
						nodes[j]->resistance = 1.; // this will be multiplied with is_resistance i.e. there is no leakage
						nodes[j]->is_resistance = 0.;
					}
					else
					{
						nodes[j]->resistance = R; // this will be multiplied with is_resistance i.e. there is leakage
						nodes[j]->is_resistance = 1.;
					}
				}
				nodes[j]->type = "junction";
				nodes[j]->type_code = 0;
				j++;
			}
			else if(sv[0] == "perif" || sv[0] == "periferia") // nodes
			{
				nodes.push_back(new moc_node(sv[1]));
				if(sv[3] == "")
				{
					nodes[j]->resistance = 1.; // this will be multiplied with is_resistance i.e. there is no leakage
					nodes[j]->is_resistance = 0.;
				}
				else
				{
					double R = stod(sv[3],0);
					if(R == 0)
					{
						nodes[j]->resistance = 1.;
						nodes[j]->is_resistance = 0.;
					}
					else
					{
						nodes[j]->resistance = R;
						nodes[j]->is_resistance = 1.;
					}
				}

				nodes[j]->pressure_out = atmospheric_pressure;
				nodes[j]->type = "periferia";
				nodes[j]->type_code = 1;
				j++;
			}
			else if(sv[0] == "perifPC" || sv[0] == "periferia_pc") // nodes
			{
				nodes.push_back(new moc_node(sv[1]));
				nodes[j]->resistance = 1.;
				nodes[j]->is_resistance = 0.;
				nodes[j]->pressure_out = stod(sv[3],0)*mmHg_to_Pa + atmospheric_pressure;
				nodes[j]->type = "periferia";
				nodes[j]->type_code = 1;
				j++;
			}
			else if(sv[0] == "sziv" || sv[0] == "heart") // nodes
			{
				// nodes, only heart
				nodes.push_back(new moc_node(sv[1]));
				nodes[j]->resistance = 1.;
				nodes[j]->is_resistance = 0.;
				nodes[j]->type = "heart";
				nodes[j]->type_code = 2;

				// loading the pressure-time curve
				if(sv.size()>4 && sv[4] != "")
				{
					string pt_file_name = sv[4];
					load_pt_series(pt_file_name);
				}

				j++;
			}
		}
	}
	else
	{
		cout << "! ERROR !" << endl << " File is not open when calling load_system_csv() function!!! file: " << file_name << "\nExiting..." << endl;
		exit(-1);
	}

	file_in.close();
}

//--------------------------------------------------------------
void solver_moc::load_pt_series(string file_name)
{
	ifstream pt_file_in;
	file_name = input_folder_path + '/' + file_name + ".csv";
	pt_file_in.open(file_name);
	if(pt_file_in.is_open())
	{
		string pt_line;
		while(getline(pt_file_in,pt_line))
		{
			vector<string> pt_sv = separate_line(pt_line);
			time_upstream.push_back(stod(pt_sv[0],0));

			double p = stod(pt_sv[1],0);
			pressure_upstream.push_back(p);
		}
		// eliminating the delay if present
		for(unsigned int i=0; i<time_upstream.size(); i++)
		{
			time_upstream[i] -= time_upstream[0];
		}
	}
	else
	{
		cout << "! ERROR !" << endl << " File is not open when calling load_pt_series() function!!! file: " << file_name << "\nExiting..." << endl;
		exit(-1);
	}
}

//--------------------------------------------------------------
void solver_moc::save_results()
{
	save_results(name);
}

//--------------------------------------------------------------
void solver_moc::save_results(string folder_name)
{
	vector<string> edge_list, node_list;
	for(int i=0; i<number_of_edges; i++)
	{
		edge_list.push_back(edges[i]->name);
	}
	for(int i=0; i<number_of_nodes; i++)
	{
		node_list.push_back(nodes[i]->name);
	}
	save_results(folder_name,edge_list,node_list);
}

//--------------------------------------------------------------
void solver_moc::save_results(string folder_name, vector<string> edge_list, vector<string> node_list)
{
   // LINUX
   mkdir("results",0777);
   mkdir(("results/" + folder_name).c_str(),0777);

   // FOR WINDOWS
   //mkdir("results");
   //mkdir(("results/" + folder_name).c_str());

   folder_name = "results/" + folder_name + "/";

   FILE *out_file;

   for(unsigned int i=0; i<node_list.size(); i++)
   {
   	int idx = node_id_to_index(node_list[i]);
      string file_name = folder_name + nodes[idx]->name + ".txt";
      out_file = fopen(file_name.c_str(),"w");
      for(unsigned int j=0; j<nodes[idx]->pressure.size(); j++)
      {
         fprintf(out_file, "%9.7e, %9.7e, %9.7e\n", time[j], nodes[idx]->pressure[j],nodes[idx]->volume_flow_rate[j]);
      }
      fclose(out_file);
   }

   for(unsigned int i=0; i<edge_list.size(); i++)
   {
   	int idx = edge_id_to_index(edge_list[i]);
      string file_name = folder_name + edges[idx]->name + ".txt";
      out_file = fopen(file_name.c_str(),"w");
      for(unsigned int j=0; j<edges[idx]->pressure_start.size(); j++)
      {
         double t = time[j];
         double ps = edges[idx]->pressure_start[j];
         double pe = edges[idx]->pressure_end[j];
         double vs = edges[idx]->velocity_start[j];
         double ve = edges[idx]->velocity_end[j];
         double vfrs = edges[idx]->volume_flow_rate_start[j];
         double vfre = edges[idx]->volume_flow_rate_end[j];
         double mfrs = edges[idx]->mass_flow_rate_start[j];
         double mfre = edges[idx]->mass_flow_rate_end[j];
         double ds = edges[idx]->diameter_start[j];
         double de = edges[idx]->diameter_end[j];
         double epszs = edges[idx]->total_deformation_start[j];
         double epsze = edges[idx]->total_deformation_end[j];
         double epsz2s = edges[idx]->damper_deformation_start[j];
         double epsz2e = edges[idx]->damper_deformation_end[j];
         double as = edges[idx]->wave_velocity_start[j];
         double ae = edges[idx]->wave_velocity_end[j];

         fprintf(out_file, "%9.7e, %9.7e, %9.7e, %9.7e, %9.7e, %9.7e, %9.7e, %9.7e, %9.7e, %9.7e, %9.7e, %9.7e, %9.7e, %9.7e, %9.7e, %9.7e, %9.7e\n",t,ps,pe,vs,ve,vfrs,vfre,mfrs,mfre,ds,de,epszs,epsze,epsz2s,epsz2e,as,ae);
      }
      fclose(out_file);
   }
}