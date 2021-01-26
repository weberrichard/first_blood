#include "../../source/moc_solver.h"

using namespace std;

int main()
{
  	srand((unsigned int) time(0));
   clock_t ido = clock();

	string case_folder = "../../models/test/";
	string case_name = "P045W_mod.csv";

  	ido = clock();
   moc_solver *patient = new moc_solver(case_folder + case_name, 5.);
	cout << endl << "\nLoading:   " << double(clock()-ido)/ CLOCKS_PER_SEC << " s" << endl;

  	ido = clock();
   patient->initialization();
	cout << endl << "\nInit:   " << double(clock()-ido)/ CLOCKS_PER_SEC << " s" << endl;

  	ido = clock();
   patient->forward_solver("C");
	cout << endl << "\nSolver:   " << double(clock()-ido)/ CLOCKS_PER_SEC << " s" << endl;

  	ido = clock();
   patient->save_results();
	cout << endl << "\nSaving:   " << double(clock()-ido)/ CLOCKS_PER_SEC << " s" << endl;

   cout << endl << endl;
   return 0;
}