#include "../../source/moc_solver.h"

using namespace std;

int main()
{
	string case_folder = "../../models/test/";
	string case_name = "P045W_mod.csv";
   moc_solver *patient = new moc_solver(case_folder + case_name, 5.);
   patient->initialization();
   patient->print_input();
   patient->forward_solver("C");
   patient->save_results();

   cout << endl << endl;
   return 0;
}