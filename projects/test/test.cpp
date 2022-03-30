#include "../../source/first_blood.h"
#include <string>

using namespace std;

vector<double> vfr_bc(vector<double> t, double q_input_f);

int main(int argc, char* argv[])
{
   // basic stuff
	string case_folder = "../../models/";
   vector<string> case_name{"Reymond_99_heart"};

   double save_dt = 1e-3;
   double heart_rate = 72.;  // if there is a heart model
   double sim_time = 20.*60./heart_rate;

   for(int i=0; i<case_name.size(); i++)
   {  
      cout << "[*] case: " << case_name[i] << endl;
      // loading original case
      first_blood *fb = new first_blood(case_folder + case_name[i]);
      fb->run();
      fb->save_results();
   }

   return 0;
}
