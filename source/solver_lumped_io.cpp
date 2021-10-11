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

			if(sv[0] == "resistor" || sv[0] == "capacitor" || sv[0] == "inductor" || sv[0] == "voltage" || sv[0] == "diode") // edges with parameter
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
				if(sv[0] == "diode")
				{
					edges[ne]->parameter = stod(sv[5],0);
					edges[ne]->type_code = 5;
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
				nodes[nn]->name = sv[1];
				nodes[nn]->pressure_initial = stod(sv[2],0);
				nodes[nn]->is_ground = false;
				nn++;
			}
			else if(sv[0] == "ground") // node with ground
			{
				nodes.push_back(new node);
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

   string fn = "results/" + folder_name + "/";

   FILE *out_file;

   for(unsigned int i=0; i<node_list.size(); i++)
   {
   	int idx = node_id_to_index(node_list[i]);
      string file_name = fn + nodes[idx]->name + ".txt";
      out_file = fopen(file_name.c_str(),"w");
      for(unsigned int j=0; j<nodes[idx]->pressure.size(); j++)
      {
      	double t = time[j];
      	double p = nodes[idx]->pressure[j];
         fprintf(out_file, "%9.7e, %9.7e\n", t, p);
      }
      fclose(out_file);
   }

   for(unsigned int i=0; i<edge_list.size(); i++)
   {
   	int idx = edge_id_to_index(edge_list[i]);
      string file_name = fn + edges[idx]->name + ".txt";
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

   string fn = "results/" + folder_name + "/";

   FILE *out_file;

   for(unsigned int i=0; i<node_list.size(); i++)
   {
   	int idx = node_id_to_index(node_list[i]);
      string file_name = fn + nodes[idx]->name + ".txt";
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

   for(unsigned int i=0; i<edge_list.size(); i++)
   {
   	int idx = edge_id_to_index(edge_list[i]);
      string file_name = fn + edges[idx]->name + ".txt";
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
