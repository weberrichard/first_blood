#include "../../source/moc_solver.h"

using namespace std;

int main()
{
	string case_folder = "../../models/test/";
	string case_name = "2cso.csv";
   moc_solver *patient = new moc_solver(case_folder + case_name);
   patient->initialization();
   patient->forward_solver("C",5.);
   patient->save_results();

   cout << endl << endl;
   return 0;
}