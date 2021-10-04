#include "../../source/moc_solver.h"

using namespace std;

int main(int argc, char* argv[])
{	
	// nodes backward calc from
	string node_1 = "cj";
	string node_2 = "rb";

	// arg inputs
	string case_name;
	double sim_time; // total time of simulation
	double alfa1, alfa2; // PWV = alfa1 + alfa2*ll_rel^2, E1 = rho*PWV*d/s
	double gamma; // E2 = gamma*E1

	if(argc == 6)
	{
		case_name = argv[1];
		sim_time = atof(argv[2]);
		alfa1 = atof(argv[3]);
		alfa2 = atof(argv[4]);
		gamma = atof(argv[5]);
	}
	else
	{
		cout << " ERROR: number of input arguments(" << argc << ") is incorrect. Proper number: 6, that is [./opt.out case_name sim_time alfa1 alfa2 gamma]" << " EXITING..." << endl;
		exit(-1);
	}

	// loading model
	//string case_folder = "../../models/";
   //moc_solver *patient = new moc_solver(case_folder + "P045.csv");
   moc_solver *patient = new moc_solver("P045.csv");

   // loading ll_rel for parameter adjustment from original case file
	vector<double> ll_rel;
	ifstream file_ll_rel;
	file_ll_rel.open(patient->input_file_path);
	string line;
	if(file_ll_rel.is_open())
	{
		while(getline(file_ll_rel,line))
		{
			vector<string> sv = patient->separate_line(line);
			
			if(sv[0] == "vis" || sv[0] == "visM") // edges
			{
				ll_rel.push_back(stod(sv[15],0));
			}
		}
	}
	file_ll_rel.close();

	// adjusting model parameters
	for(unsigned int i=0; i<patient->edges.size(); i++)
	{
		double a = alfa1 + alfa2*ll_rel[i]*ll_rel[i];
		double dn = patient->edges[i]->nominal_diameter;
		double sn = patient->edges[i]->nominal_thickness;
		double rho = patient->density;
		double E1 =  rho * a*a * dn / sn;
		patient->edges[i]->elasticity_spring = E1;
		double E2 = gamma * E1;
		patient->edges[i]->elasticity_voigt = E2;
	}

	// ----------------------- //
	// full solver from node_1 //
	// ----------------------- //

	// modifying input time-series
	patient->time_upstream.clear();
	patient->pressure_upstream.clear();
	ifstream file_pres;
	file_pres.open("measurements/" + case_name + "/" + case_name + "_" + node_1 + ".csv");
	if(file_pres.is_open())
	{
		while(getline(file_pres,line))
		{
			vector<string> sv = patient->separate_line(line);
			
			double t = stod(sv[0],0);
			patient->time_upstream.push_back(t);

			double p = stod(sv[1],0);
			p *= patient->mmHg_to_Pa; // converting the units from mmHg to Pa
			p += patient->atmospheric_pressure; // adding the atmospheric pressure
			patient->pressure_upstream.push_back(p);
		}
	}
	file_pres.close();

	vector<double> tu = patient->time_upstream;
	vector<double> pu = patient->pressure_upstream;

   patient->initialization();
   patient->full_solver(node_1,sim_time);
   vector<string> node_list, edge_list;
   node_list.push_back("C");
   patient->save_results(case_name + "_" + node_1, edge_list, node_list);

	// ----------------------- // 
	// full solver from node_2 //
	// ----------------------- // 

   patient->time_upstream = tu;
   patient->pressure_upstream = pu;

	// modifying input time-series
	patient->time_upstream.clear();
	patient->pressure_upstream.clear();
	file_pres.open("measurements/" + case_name + "/" + case_name + "_" + node_2 + ".csv");
	if(file_pres.is_open())
	{
		while(getline(file_pres,line))
		{
			vector<string> sv = patient->separate_line(line);
			
			patient->time_upstream.push_back(stod(sv[0],0));

			double p = stod(sv[1],0);
			p *= patient->mmHg_to_Pa; // converting the units from mmHg to Pa
			p += patient->atmospheric_pressure; // adding the atmospheric pressure
			patient->pressure_upstream.push_back(p);
		}
	}
	file_pres.close();

   patient->initialization();
   patient->full_solver(node_2,sim_time);
   patient->save_results(case_name + "_" + node_2, edge_list, node_list);

   cout << endl << endl;
   return 0;
}