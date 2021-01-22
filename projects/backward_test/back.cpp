#include "../../source/moc_solver.h"

using namespace std;

int main()
{
	srand((unsigned int) time(0));
   clock_t ido = clock();

	string case_folder = "../../models/test/";
	string case_name = "P045.csv";
   moc_solver *patient = new moc_solver(case_folder + case_name);

  	ido = clock();
   patient->forward_solver("C",5.);
	cout << endl << "\n Forward solver:   " << double(clock()-ido)/ CLOCKS_PER_SEC << " s" << endl;
   patient->save_results();

   // saving inner node time, pressure
	ofstream out_file;
	out_file.open (case_folder + "ger_back.csv");
	int edge_idx = patient->edge_id_to_index("A42");
	for(unsigned int i=0; i<patient->time.size(); i++)
	{
		out_file << patient->time[i] << ", " << (patient->edges[edge_idx]->pressure_start[i]-patient->atmospheric_pressure)/patient->mmHg_to_Pa << "\n";
	}
	out_file.close();

	cout << endl << endl;

	case_name = "P045_back.csv";
   moc_solver *patient_back = new moc_solver(case_folder + case_name);

  	ido = clock();
   patient_back->full_solver("21a",5.);
	cout << endl << "\n Full solver:   " << double(clock()-ido)/ CLOCKS_PER_SEC << " s" << endl;
   patient_back->save_results();

   cout << endl << endl;
   return 0;
}