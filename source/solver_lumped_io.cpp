#include "solver_lumped.h"

//--------------------------------------------------------------
void solver_lumped::load_model()
{
	ifstream file_in;
	string file_name = input_folder_path + '/' + name + ".csv";
	file_in.open(file_name);
	string line;
	if(file_in.is_open())
	{
		int nn=0,ne=0; // ne for edges, nn for nodes
		while(getline(file_in,line))
		{	
			// cleaning unnecessary characters
			line.erase(remove(line.begin(), line.end(), ' '), line.end());
			line.erase(remove(line.begin(), line.end(), '\n'), line.end());
			line.erase(remove(line.begin(), line.end(), '\r'), line.end());
			
			vector<string> sv = separate_line(line);

			if(sv[0] == "resistor" || sv[0] == "capacitor" || sv[0] == "inductor" || sv[0] == "voltage" || sv[0] == "diode" || sv[0] == "resistor2" || sv[0] == "valve" || sv[0] == "resistor_coronary" || sv[0] == "capacitor_coronary" || sv[0] == "current") // edges with parameter
			{
				edges.push_back(new edge);
				edges[ne]->type = sv[0];
				edges[ne]->name = sv[1];
				edges[ne]->node_name_start = sv[2];
				edges[ne]->node_name_end = sv[3];
				edges[ne]->volume_flow_rate_initial = stod(sv[4],0);
				if(sv[0] == "resistor")
				{
					edges[ne]->parameter = stod(sv[5],0);
					edges[ne]->type_code = 0;
				}
				else if(sv[0] == "capacitor")
				{
					edges[ne]->parameter = stod(sv[5],0);
					edges[ne]->type_code = 1;
				}
				else if(sv[0] == "inductor")
				{
					edges[ne]->parameter = stod(sv[5],0);
					edges[ne]->type_code = 3;
				}
				else if(sv[0] == "voltage")
				{
					edges[ne]->parameter = stod(sv[5],0);					
					edges[ne]->type_code = 4;
				}
				else if(sv[0] == "diode")
				{
					edges[ne]->parameter = stod(sv[5],0);
					edges[ne]->type_code = 5;
				}
				else if(sv[0] == "resistor2" || sv[0] == "valve")
				{
					edges[ne]->parameter = stod(sv[5],0);
					edges[ne]->type_code = 6;
				}
				else if(sv[0] == "resistor_coronary")
				{
					edges[ne]->parameter = stod(sv[5],0);
					edges[ne]->type_code = 7;
				}
				else if(sv[0] == "capacitor_coronary")
				{
					edges[ne]->parameter = stod(sv[5],0);
					edges[ne]->type_code = 8;
				}
				else if(sv[0] == "current")
				{
					edges[ne]->parameter = stod(sv[5],0);					
					edges[ne]->type_code = 9;
				}
				ne++;
			}
			else if(sv[0] == "elastance") // edges without parameters
			{
				edges.push_back(new edge);
				edges[ne]->type = sv[0];
				edges[ne]->name = sv[1];
				edges[ne]->node_name_start = sv[2];
				edges[ne]->node_name_end = sv[3];
				edges[ne]->volume_flow_rate_initial = stod(sv[4],0);
				edges[ne]->type_code = 2;
				ne++;
			}
			else if(sv[0] == "node") // node
			{
				nodes.push_back(new node);
				nodes[nn]->type = sv[0];
				nodes[nn]->name = sv[1];
				nodes[nn]->pressure_initial = stod(sv[2],0);
				nodes[nn]->is_ground = false;
				nn++;
			}
			else if(sv[0] == "ground") // node with ground
			{
				nodes.push_back(new node);
				nodes[nn]->type = sv[0];
				nodes[nn]->name = sv[1];
				nodes[nn]->pressure_initial = stod(sv[2],0);
				nodes[nn]->is_ground = true;
				nn++;
			}
		}
	}
	else
	{
		std::cout << "! ERROR !" << endl << " File is not open when calling load_system_csv() function!!! file: " << file_name << "\nExiting..." << endl;
		exit(-1);
	}

	// setting size of elements
	number_of_nodes = nodes.size();
	number_of_edges = edges.size();

	file_in.close();
}

//--------------------------------------------------------------
void solver_lumped::save_results()
{
	save_results(name);
}

//--------------------------------------------------------------
void solver_lumped::save_results(string folder_name)
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
void solver_lumped::save_results(string folder_name, vector<string> edge_list, vector<string> node_list)
{
   // LINUX
   mkdir("results",0777);
   mkdir(("results/" + folder_name).c_str(),0777);

   string fn = "results/" + folder_name + "/" + name;

   FILE *out_file;

   for(unsigned int i=0; i<node_list.size(); i++)
   {
   	int idx = node_id_to_index(node_list[i]);
   	if(nodes[idx]->do_save_memory)
   	{		
	   	mkdir(fn.c_str(),0777);
		   string file_name = fn + "/" + nodes[idx]->name + ".txt";
	      out_file = fopen(file_name.c_str(),"w");
	      for(unsigned int j=0; j<nodes[idx]->pressure.size(); j++)
	      {
	      	double t = time[j];
	      	double p = nodes[idx]->pressure[j];
	         fprintf(out_file, "%9.7e, %9.7e\n", t, p);
	      }
	      fclose(out_file);
   	}
   }

   for(unsigned int i=0; i<edge_list.size(); i++)
   {
   	int idx = edge_id_to_index(edge_list[i]);
   	if(edges[idx]->do_save_memory)
   	{
	   	mkdir(fn.c_str(),0777);
		   string file_name = fn + "/" + edges[idx]->name + ".txt";
		   out_file = fopen(file_name.c_str(),"w");
		   for(unsigned int j=0; j<edges[idx]->volume_flow_rate.size(); j++)
		   {
		      double t = time[j];
		      double vfr = edges[idx]->volume_flow_rate[j];
		      fprintf(out_file, "%9.7e, %9.7e\n",t,vfr);
		   }
		   fclose(out_file);
   	}
   }
}

//--------------------------------------------------------------
void solver_lumped::save_results(double dt, string folder_name)
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
	save_results(dt,folder_name,edge_list,node_list);
}

//--------------------------------------------------------------
void solver_lumped::save_results(double dt, string folder_name, vector<string> edge_list, vector<string> node_list)
{
   // LINUX
   mkdir("results",0777);
   mkdir(("results/" + folder_name).c_str(),0777);
   
   string fn = "results/" + folder_name + "/" + name;

   FILE *out_file;

   for(unsigned int i=0; i<node_list.size(); i++)
   {
   	int idx = node_id_to_index(node_list[i]);
   	if(nodes[idx]->do_save_memory)
   	{
   		mkdir(fn.c_str(),0777);
	      string file_name = fn + "/" + nodes[idx]->name + ".txt";
	      out_file = fopen(file_name.c_str(),"w");

	      int j=0;
			double ts=0.;
			double t_end=time.back();
			while(ts<t_end && j<time.size()-1)
			{
				if(time[j]<=ts && ts<time[j+1])
				{	
					double a0 = (time[j+1]-ts)/(time[j+1]-time[j]);
					double a1 = (ts-time[j])/(time[j+1]-time[j]);

			   	double p = nodes[idx]->pressure[j]*a0 + nodes[idx]->pressure[j+1]*a1;
			      fprintf(out_file, "%9.7e, %9.7e\n", ts, p);
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

   for(unsigned int i=0; i<edge_list.size(); i++)
   {
   	int idx = edge_id_to_index(edge_list[i]);
   	if(edges[idx]->do_save_memory)
   	{
   		mkdir(fn.c_str(),0777);
	      string file_name = fn + "/" + edges[idx]->name + ".txt";
	      out_file = fopen(file_name.c_str(),"w");

	      int j=0;
			double ts=0.;
			double t_end=time.back();
			while(ts<t_end && j<time.size()-1)
			{
				if(time[j]<=ts && ts<time[j+1])
				{
		         double a0 = (time[j+1]-ts)/(time[j+1]-time[j]);
					double a1 = (ts-time[j])/(time[j+1]-time[j]);

		         double vfr = edges[idx]->volume_flow_rate[j]*a0 + edges[idx]->volume_flow_rate[j+1]*a1;
		         fprintf(out_file, "%9.7e, %9.7e\n",ts,vfr);
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

//--------------------------------------------------------------
void solver_lumped::save_model(string model_name, string folder_name)
{
   FILE *out_file;
   string file_name = folder_name + "/" + model_name + "/" + name + ".csv";
	out_file = fopen(file_name.c_str(),"w");

	// writing edges data to file
	fprintf(out_file,"data of edges\n");
	fprintf(out_file,"type, name, node start, node end, initial condition [SI], parameter [SI]\n");
	for(int i=0; i<number_of_edges; i++)
	{
		fprintf(out_file, "%s, %s, %s, %s, %6.3e, %6.3e\n",edges[i]->type.c_str(),edges[i]->name.c_str(),edges[i]->node_name_start.c_str(), edges[i]->node_name_end.c_str(), edges[i]->volume_flow_rate_initial, edges[i]->parameter);
	}
	fprintf(out_file,"\n");

	// writing nodes data to file
	fprintf(out_file,"data of nodes\n");
	fprintf(out_file,"type, name, initial condition [SI]\n");
	for(int i=0; i<number_of_nodes; i++)
	{
		fprintf(out_file, "%s, %s, %6.3e\n",nodes[i]->type.c_str(),nodes[i]->name.c_str(),nodes[i]->pressure_initial);
	}
	fprintf(out_file,"\n");
   fclose(out_file);
}

//--------------------------------------------------------------
void solver_lumped::save_initials(string model_name, string folder_name)
{
   string file_name = folder_name + "/" + model_name + "/init/" + name + ".csv";

	FILE* out_file;
	out_file = fopen(file_name.c_str(),"w");

	for(int i=0; i<number_of_edges; i++)
	{	
		double vfr = edges[i]->vfr/1.e6; // convert ml/s to m3/s
		fprintf(out_file, "edge, %s, %9.7e\n",edges[i]->name.c_str(),vfr);
	}

	for(int i=0; i<number_of_nodes; i++)
	{
		double p = nodes[i]->p*mmHg_to_Pa; // convert mmHg to Pa
		fprintf(out_file, "node, %s, %9.7e\n",nodes[i]->name.c_str(),p);
	}

   fclose(out_file);
}

//--------------------------------------------------------------
void solver_lumped::load_initials()
{
	ifstream file_in;
	string file_name = input_folder_path + "/init/" + name + ".csv";
	file_in.open(file_name);
	string line;
	if(file_in.is_open())
	{
		double p_conv = 1./(mmHg_to_Pa*elastance(0.));
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
				if(sv[0] == "edge")
				{
					int idx = edge_id_to_index(sv[1]);
					double vfr = stod(sv[2],0); // m3/s
					edges[idx]->volume_flow_rate.clear();
					edges[idx]->volume_flow_rate.push_back(vfr);
					edges[idx]->vfr = vfr*1.e6; // to ml/s
				}
				else if(sv[0] == "node")
				{
					int idx = node_id_to_index(sv[1]);
					double p = stod(sv[2],0); // Pa
					nodes[idx]->pressure.clear();
					nodes[idx]->pressure.push_back(p);
					nodes[idx]->p = p/mmHg_to_Pa; // to mmHg
					nodes[idx]->y = p*p_conv;
				}
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
