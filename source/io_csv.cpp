#include "first_blood.h"

using namespace std;

//--------------------------------------------------------------
void first_blood::load_system_csv()
{
	ifstream file_in;
	file_in.open(input_file_path);
	string line;
	if(file_in.is_open())
	{
		int i=0,j=0; // i for edges, j for nodes
		while(getline(file_in,line))
		{
			vector<string> sv = separate_line(line);
			
			if(sv[0] == "vis" || sv[0] == "visM") // edges
			{
				edges.push_back(new edge(sv[1]));
				edges[i]->node_name_start = sv[2];
				edges[i]->node_name_end = sv[3];
				edges[i]->nominal_diameter = stod(sv[4],0);
				edges[i]->nominal_thickness = stod(sv[5],0);
				edges[i]->length = stod(sv[6],0);
				edges[i]->division_points = stoi(sv[7],0);
				edges[i]->elasticity_spring = stod(sv[8],0);
				edges[i]->elasticity_voigt = stod(sv[9],0);
				edges[i]->viscosity = stod(sv[10],0);
				edges[i]->resistance_start = stod(sv[11],0);
				edges[i]->resistance_end = stod(sv[12],0);
				i++;
			}
			else if(sv[0] == "elag" || sv[0] == "junction") // nodes
			{
				nodes.push_back(new node(sv[1]));
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
						nodes[j]->resistance = R; // this will be multiplied with is_resistance i.e. there is no leakage
						nodes[j]->is_resistance = 1.;
					}
				}
				nodes[j]->pressure_out = atmospheric_pressure;
				nodes[j]->type = "junction";
				nodes[j]->type_code = 0;
				j++;
			}
			else if(sv[0] == "perif" || sv[0] == "periferia") // nodes
			{
				nodes.push_back(new node(sv[1]));
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
						nodes[j]->resistance = R; // this will be multiplied with is_resistance i.e. there is no leakage
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
				nodes.push_back(new node(sv[1]));
				nodes[j]->resistance = 1.;
				nodes[j]->is_resistance = 0.;
				nodes[j]->pressure_out = stod(sv[3],0)*mmHg_to_Pa + atmospheric_pressure;
				nodes[j]->type = "periferia";
				nodes[j]->type_code = 1;
				j++;
			}
			else if(sv[0] == "sziv") // nodes
			{
				// nodes, only heart
				nodes.push_back(new node(sv[1]));
				nodes[j]->resistance = 1.; // this will be multiplied with is_resistance i.e. there is no leakage
				nodes[j]->is_resistance = 0.;
				nodes[j]->is_upstream_boundary = true;
				nodes[j]->type = "heart";
				nodes[j]->type_code = 2;

				// loading the pressure-time curve
				string pt_file_name = sv[4] + ".csv";
				ifstream pt_file_in;
				pt_file_in.open(input_folder_path + pt_file_name);
				if(pt_file_in.is_open())
				{
					string pt_line;
					while(getline(pt_file_in,pt_line))
					{
						vector<string> pt_sv = separate_line(pt_line);
						time_upstream.push_back(stod(pt_sv[0],0));

						double p = stod(pt_sv[1],0);
         			p *= mmHg_to_Pa; // converting the units from mmHg to Pa
         			p += atmospheric_pressure; // adding the atmospheric pressure

						pressure_upstream.push_back(p);
					}
				}
				else
				{
					cout << "! ERROR !" << endl << " File is not open when calling load_system_csv() function!!! file: " << input_folder_path + pt_file_name << "\nExiting..." << endl;
					exit(-1);
				}

				j++;
			}
		}
	}
	else
	{
		cout << "! ERROR !" << endl << " File is not open when calling load_system_csv() function!!! file: " << input_file_path << "\nExiting..." << endl;
		exit(-1);
	}

	file_in.close();
}

//--------------------------------------------------
vector<string> first_blood::separate_line(string line)
{
	string s="";
	vector<string> sv;
	for(string::iterator i=line.begin(); i!=line.end(); i++)
	{
		if(*i!=',')
		{
			s += *i;
		}
		else
		{
			sv.push_back(s);
			s="";
		}
	}
	sv.push_back(s);

	return sv;
}


//--------------------------------------------------------------
void first_blood::save_results()
{
   // LINUX
   mkdir("results",0777);
   mkdir(("results/" + case_name).c_str(),0777);

   // FOR WINDOWS
   //mkdir(case_name.c_str());

   string folder_name = "results/" + case_name + "/";

   FILE *out_file;
   for(unsigned int i=0; i<number_of_nodes; i++)
   {  
      string file_name = folder_name + nodes[i]->name + ".txt";
      out_file = fopen(file_name.c_str(),"w");
      for(unsigned int j=0; j<number_of_timesteps; j++)
      {
         fprintf(out_file, "%9.7e, %9.7e, %9.7e\n", time[j], nodes[i]->pressure[j],nodes[i]->volume_flow_rate[j]);
      }
      fclose(out_file);
   }
   for(unsigned int i=0; i<number_of_edges; i++)
   {
      string file_name = folder_name + edges[i]->name + ".txt";
      out_file = fopen(file_name.c_str(),"w");
      for(unsigned int j=0; j<number_of_timesteps; j++)
      {
         double t = time[j];
         double ps = edges[i]->pressure_start[j];
         double pe = edges[i]->pressure_end[j];
         double vs = edges[i]->velocity_start[j];
         double ve = edges[i]->velocity_end[j];
         double vfrs = edges[i]->volume_flow_rate_start[j];
         double vfre = edges[i]->volume_flow_rate_end[j];
         double mfrs = edges[i]->mass_flow_rate_start[j];
         double mfre = edges[i]->mass_flow_rate_end[j];
         double ds = edges[i]->diameter_start[j];
         double de = edges[i]->diameter_end[j];
         double epszs = edges[i]->total_deformation_start[j];
         double epsze = edges[i]->total_deformation_end[j];
         double epsz2s = edges[i]->damper_deformation_start[j];
         double epsz2e = edges[i]->damper_deformation_end[j];
         double as = edges[i]->wave_velocity_start[j];
         double ae = edges[i]->wave_velocity_end[j];

         fprintf(out_file, "%9.7e %9.7e %9.7e %9.7e %9.7e %9.7e %9.7e %9.7e %9.7e %9.7e %9.7e %9.7e %9.7e %9.7e %9.7e %9.7e %9.7e\n",t,ps,pe,vs,ve,vfrs,vfre,mfrs,mfre,ds,de,epszs,epsze,epsz2s,epsz2e,as,ae);
      }
      fclose(out_file);
   }
}