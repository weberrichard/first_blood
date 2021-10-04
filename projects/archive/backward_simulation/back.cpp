#include "../../source/moc_solver.h"

using namespace std;

int main(int argc, char* argv[])
{
	string case_name;
	double sim_time;
	string input_node;
	if(argc == 1)
	{
		case_name =  "4_pipe.csv";
		sim_time = 5.;
		input_node = "N2";
	}
	else if(argc == 2)
	{
		case_name = argv[1];
		sim_time = 5.;
		input_node = "N2";
	}
	else if(argc == 3)
	{
		case_name = argv[1];
		sim_time = atof(argv[2]);
		input_node = "N2";
	}
	else if(argc == 4)
	{
		case_name = argv[1];
		sim_time = atof(argv[2]);
		input_node = argv[3];
	}

	srand((unsigned int) time(0));
   clock_t ido = clock();

	string case_folder = "../../models/";
   moc_solver *patient_back = new moc_solver(case_folder + case_name);

  	ido = clock();
   patient_back->full_solver(input_node,sim_time);
	cout << endl << "\n Full solver:   " << double(clock()-ido)/ CLOCKS_PER_SEC << " s" << endl;
   patient_back->save_results();

   cout << endl << endl;
   return 0;
}