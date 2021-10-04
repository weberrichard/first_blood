#include "../../source/moc_solver.h"

using namespace std;

int main(int argc, char* argv[])
{
	string case_name;
	double sim_time;
	if(argc == 1)
	{
		case_name =  "4_pipe.csv";
		sim_time = 5.;
	}
	else if(argc == 2)
	{
		case_name = argv[1];
		sim_time = 5.;
	}
	else if(argc == 3)
	{
		case_name = argv[1];
		sim_time = atof(argv[2]);
	}

	string case_folder = "../../models/";
   moc_solver *patient = new moc_solver(case_folder + case_name);
   patient->initialization();
   patient->forward_solver("C",sim_time);
   patient->save_results();

   cout << endl << endl;
   return 0;
}