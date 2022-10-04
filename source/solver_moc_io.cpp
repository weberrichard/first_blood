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
				edges[i]->type = sv[0];
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
				nodes[j]->type = "node";
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
				nodes[j]->type = "perif";
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

				if(sv.size()>3)
				{
					if(sv[3] == "P")
					{
						type_upstream = 0;
					}
					else if(sv[3] == "Q")
					{
						type_upstream = 1;
					}
					else if(sv[3] == "V")
					{
						type_upstream = 2;
					}
				}
				// loading the pressure-time curve
				if(sv.size()>4 && sv[4] != "")
				{
					pt_file_name = sv[4];
					load_time_series(pt_file_name);
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

	// setting size of elements
	number_of_nodes = nodes.size();
	number_of_edges = edges.size();

	file_in.close();
}

//--------------------------------------------------------------
void solver_moc::load_time_series(string file_name)
{
	time_upstream.clear();
	value_upstream.clear();
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
			value_upstream.push_back(p);
		}
		// eliminating the delay if present
		for(unsigned int i=0; i<time_upstream.size(); i++)
		{
			time_upstream[i] -= time_upstream[0];
		}
	}
	else
	{
		cout << "! ERROR !" << endl << " File is not open when calling load_time_series() function!!! file: " << file_name << "\nExiting..." << endl;
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
		edge_list.push_back(edges[i]->ID);
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
   
	folder_name = folder_name + "/" + name;
   mkdir(("results/" + folder_name).c_str(),0777);
   
   folder_name = "results/" + folder_name + "/";

   FILE *out_file;

   for(unsigned int i=0; i<node_list.size(); i++)
   {
   	int idx = node_id_to_index(node_list[i]);
   	if(idx>-1)
   	{
   		if(nodes[idx]->do_save_memory)
			{
		      string file_name = folder_name + nodes[idx]->name + ".txt";
		      out_file = fopen(file_name.c_str(),"w");
		      for(unsigned int j=0; j<nodes[idx]->time.size(); j++)
		      {
		         fprintf(out_file, "%9.7e, %9.7e, %9.7e\n", nodes[idx]->time[j], nodes[idx]->pressure[j],nodes[idx]->volume_flow_rate[j]);
		      }
		      fclose(out_file);
			}
   	}
   }

   for(unsigned int i=0; i<edge_list.size(); i++)
   {
   	int idx = edge_id_to_index(edge_list[i]);
   	if(idx>-1)
   	{
   		if(edges[idx]->do_save_memory)
   		{
	      	string file_name = folder_name + edges[idx]->ID + ".txt";
		      out_file = fopen(file_name.c_str(),"w");\
		      for(unsigned int j=0; j<edges[idx]->pressure_start.size(); j++)
		      {
		         double t = edges[idx]->time[j];
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
   }
}

//--------------------------------------------------------------
void solver_moc::save_results(double dt, string folder_name)
{
	vector<string> edge_list, node_list;
	for(int i=0; i<number_of_edges; i++)
	{
		edge_list.push_back(edges[i]->ID);
	}
	for(int i=0; i<number_of_nodes; i++)
	{
		node_list.push_back(nodes[i]->name);
	}
	save_results(dt,folder_name,edge_list,node_list);
}

//--------------------------------------------------------------
void solver_moc::save_results(double dt, string folder_name, vector<string> edge_list, vector<string> node_list)
{
   // LINUX
   mkdir("results",0777);
   mkdir(("results/" + folder_name).c_str(),0777);
   
	folder_name = folder_name + "/" + name;
   mkdir(("results/" + folder_name).c_str(),0777);
   folder_name = "results/" + folder_name + "/";

   FILE *out_file;

   for(unsigned int i=0; i<node_list.size(); i++)
   {
   	int idx = node_id_to_index(node_list[i]);
   	if(idx>-1)
   	{
	   	if(nodes[idx]->do_save_memory)
	   	{
		      string file_name = folder_name + nodes[idx]->name + ".txt";
		      out_file = fopen(file_name.c_str(),"w");
				
				int j=0;
				double ts=0.;
				double t_end=nodes[idx]->time.back();
				while(ts<t_end && j<nodes[idx]->time.size()-1)
				{
					if(nodes[idx]->time[j]<=ts && ts<nodes[idx]->time[j+1])
					{	
						double a0 = (nodes[idx]->time[j+1]-ts)/(nodes[idx]->time[j+1]-nodes[idx]->time[j]);
						double a1 = (ts-nodes[idx]->time[j])/(nodes[idx]->time[j+1]-nodes[idx]->time[j]);
						double pp = nodes[idx]->pressure[j]*a0 + nodes[idx]->pressure[j+1]*a1;
						double qp = nodes[idx]->volume_flow_rate[j]*a0 + nodes[idx]->volume_flow_rate[j+1]*a1;
		         	fprintf(out_file, "%9.7e, %9.7e, %9.7e\n", ts, pp, qp);
		         	ts += dt;
					}
					else
					{
						j++;
					}
				}
		      fclose(out_file);
	   	}
   	}
   }

   for(unsigned int i=0; i<edge_list.size(); i++)
   {
   	int idx = edge_id_to_index(edge_list[i]);
   	if(idx>-1)
   	{
	   	if(edges[idx]->do_save_memory)
	   	{
		      string file_name = folder_name + edges[idx]->ID + ".txt";
		      out_file = fopen(file_name.c_str(),"w");

		   	int j=0;
		   	double ts=0.;
		   	double t_end=edges[idx]->time.back();
		   	while(ts<t_end && j<edges[idx]->time.size()-1)
		   	{
					if(edges[idx]->time[j]<=ts && ts<edges[idx]->time[j+1])
					{
						double a0 = (edges[idx]->time[j+1]-ts)/(edges[idx]->time[j+1]-edges[idx]->time[j]);
						double a1 = (ts-edges[idx]->time[j])/(edges[idx]->time[j+1]-edges[idx]->time[j]);

			         double ps = edges[idx]->pressure_start[j]*a0 + edges[idx]->pressure_start[j+1]*a1;
			         double pe = edges[idx]->pressure_end[j]*a0 + edges[idx]->pressure_end[j+1]*a1;
			         double vs = edges[idx]->velocity_start[j]*a0 + edges[idx]->velocity_start[j+1]*a1;
			         double ve = edges[idx]->velocity_end[j]*a0 + edges[idx]->velocity_end[j+1]*a1;
			         double vfrs = edges[idx]->volume_flow_rate_start[j]*a0 + edges[idx]->volume_flow_rate_start[j+1]*a1;
			         double vfre = edges[idx]->volume_flow_rate_end[j]*a0 + edges[idx]->volume_flow_rate_end[j+1]*a1;
			         double mfrs = edges[idx]->mass_flow_rate_start[j]*a0 + edges[idx]->mass_flow_rate_start[j+1]*a1;
			         double mfre = edges[idx]->mass_flow_rate_end[j]*a0 + edges[idx]->mass_flow_rate_end[j+1]*a1;
			         double ds = edges[idx]->diameter_start[j]*a0 + edges[idx]->diameter_start[j+1]*a1;
			         double de = edges[idx]->diameter_end[j]*a0 + edges[idx]->diameter_end[j+1]*a1;
			         double epszs = edges[idx]->total_deformation_start[j]*a0 + edges[idx]->total_deformation_start[j+1]*a1;
			         double epsze = edges[idx]->total_deformation_end[j]*a0 + edges[idx]->total_deformation_end[j+1]*a1;
			         double epsz2s = edges[idx]->damper_deformation_start[j]*a0 + edges[idx]->damper_deformation_start[j+1]*a1;
			         double epsz2e = edges[idx]->damper_deformation_end[j]*a0 + edges[idx]->damper_deformation_end[j+1]*a1;
			         double as = edges[idx]->wave_velocity_start[j]*a0 + edges[idx]->wave_velocity_start[j+1]*a1;
			         double ae = edges[idx]->wave_velocity_end[j]*a0 + edges[idx]->wave_velocity_end[j+1]*a1;

		         	fprintf(out_file, "%9.7e, %9.7e, %9.7e, %9.7e, %9.7e, %9.7e, %9.7e, %9.7e, %9.7e, %9.7e, %9.7e, %9.7e, %9.7e, %9.7e, %9.7e, %9.7e, %9.7e\n",ts,ps,pe,vs,ve,vfrs,vfre,mfrs,mfre,ds,de,epszs,epsze,epsz2s,epsz2e,as,ae);
		         	ts += dt;
					}
					else
					{
						j++;
					}
		   	}
		      fclose(out_file);
	   	}
   	}
   }
}

//--------------------------------------------------------------
void solver_moc::save_model(string model_name, string folder_name)
{
	// saving the input pt series first
	if(time_upstream.size()>0)
	{
		save_pt_series(model_name, folder_name);
	}

	// saving the model
	FILE* out_file;
   string file_name = folder_name + "/" + model_name + "/" + name + ".csv";
	out_file = fopen(file_name.c_str(),"w");

	// saving edges
	fprintf(out_file,"type,ID,name,start_node,end_node,start_diameter[SI],end_diameter[SI],start_thickness[SI],end_thickness[SI],length[SI],division_points,elastance_1[SI],elastance_2[SI],eta[SI],res_start[SI],res_end[SI]\n");
	for(int i=0; i<number_of_edges; i++)
	{
		fprintf(out_file,"%s,%s,%s,%s,%s,%6.3e,%6.3e,%6.3e,%6.3e,%6.3e,%i,%6.3e,%6.3e,%6.3e,%6.3e,%6.3e\n",edges[i]->type.c_str(),edges[i]->ID.c_str(),edges[i]->name.c_str(),edges[i]->node_name_start.c_str(),edges[i]->node_name_end.c_str(),edges[i]->nominal_diameter_start,edges[i]->nominal_diameter_end,edges[i]->nominal_thickness_start,edges[i]->nominal_thickness_end,edges[i]->length,edges[i]->division_points,edges[i]->elasticity_spring,edges[i]->elasticity_voigt,edges[i]->viscosity,edges[i]->resistance_start,edges[i]->resistance_end);
	}
	fprintf(out_file,"\n");

	// saving nodes
	fprintf(out_file, "type,ID,name,valami,parameter,file name\n");
	for(int i=0; i<number_of_nodes; i++)
	{
		if(nodes[i]->type_code == 2) // heart
		{
			fprintf(out_file, "%s,%s,0,P,%s\n",nodes[i]->type.c_str(),nodes[i]->name.c_str(),pt_file_name.c_str());
		}
		else if(nodes[i]->type_code == 0) // node
		{
			fprintf(out_file, "%s,%s,0,%6.3e\n",nodes[i]->type.c_str(),nodes[i]->name.c_str(),nodes[i]->resistance);
		}
		else if(nodes[i]->type_code == 1) // perif
		{
			fprintf(out_file, "%s,%s,0\n",nodes[i]->type.c_str(),nodes[i]->name.c_str());
		}
	}
	fprintf(out_file,"\n");
   fclose(out_file);
}

//--------------------------------------------------------------
void solver_moc::save_pt_series(string model_name, string folder_name)
{
	FILE* out_file;
   string file_name = folder_name + "/" + model_name + "/" + pt_file_name + ".csv";
	out_file = fopen(file_name.c_str(),"w");

	if(type_upstream == 0)
	{
		for(int i=0; i<time_upstream.size(); i++)
	   {
			double p_mmHg = (value_upstream[i]-atmospheric_pressure)/mmHg_to_Pa;   	
			fprintf(out_file,"%6.3e,%8.5e\n",time_upstream[i],p_mmHg);	
	   }
	}
	else if(type_upstream == 1)
	{
		for(int i=0; i<time_upstream.size(); i++)
		{
			double q_mls = value_upstream[i]*1.e6; // m3/s to ml/s
			fprintf(out_file,"%6.3e,%8.5e\n",time_upstream[i],q_mls);	
		}
	}
   
   fclose(out_file);

}

//--------------------------------------------------------------
void solver_moc::save_initials(string model_name, string folder_name)
{
   string file_name = folder_name + "/" + model_name + "/init/" + name + ".csv";

	FILE* out_file;
	out_file = fopen(file_name.c_str(),"w");

	for(int i=0; i<number_of_edges; i++)
	{
		edges[i]->save_initials(out_file);
	}

   fclose(out_file);
}

//--------------------------------------------------------------
void solver_moc::load_initials()
{
	ifstream file_in;
	string file_name = input_folder_path + "/init/" + name + ".csv";
	file_in.open(file_name);
	string line;
	if(file_in.is_open())
	{
		while(getline(file_in,line))
		{
			// clearing spaces and \n
			line.erase(remove(line.begin(), line.end(), ' '), line.end());
			line.erase(remove(line.begin(), line.end(), '\n'), line.end());
			line.erase(remove(line.begin(), line.end(), '\r'), line.end());

			// seperating the strings by comma
			vector<string> sv = separate_line(line);

			if(sv.size()>0)
			{
				string id = sv[0];
				int idx = edge_id_to_index(id);
				vector<double> ic(sv.size()-1);
				for(int i=0; i<ic.size(); i++)
				{
					ic[i] = stod(sv[i+1],0);
				}
				edges[idx]->set_initials(ic);
			}
		}
	}
	else
	{
		cout << "! ERROR !" << endl << " File is not open when calling load_initials() function!!! file: " << file_name << "\nExiting..." << endl;
		exit(-1);
	}
	file_in.close();
}
