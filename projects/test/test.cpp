#include "../../source/first_blood.h"
#include <string>

using namespace std;

vector<double> vfr_bc(vector<double> t, double q_input_f);

int main(int argc, char* argv[])
{
   // basic stuff
   double dt=1.e-3;
	string case_folder = "../../models/";
   vector<string> case_name;
   case_name.push_back("moc_test");
   case_name.push_back("nonlinear_test");
   case_name.push_back("Reymond_99_heart");
   case_name.push_back("Reymond_99_heart_ref");
   srand((unsigned int) time(0));
   clock_t ido = clock();

   for(int i=0; i<case_name.size(); i++)
   {  
      cout << "[*] case: " << case_name[i] << endl;
      // loading original case
      first_blood *fb = new first_blood(case_folder + case_name[i]);
      cout << "LOAD OK " << endl;
      ido = clock();
      fb->run();
      cout << endl << case_name[i] << ": " << double(clock()-ido)/ CLOCKS_PER_SEC << " s" << endl;
      cout << "RUN OK " << endl;
      fb->save_results(dt);
      cout << "SAVE OK " << endl;
   }

   return 0;
}
