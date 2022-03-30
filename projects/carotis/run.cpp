#include "../../source/first_blood.h"
#include <string>

using namespace std;

vector<double> vfr_bc(vector<double> t, double q_input_f);

int main(int argc, char* argv[])
{
   // basic stuff
	string case_folder = "models/";
   string case_name = "Carotis_1";
   if(argc == 2)
   {
      case_name = argv[1];
   }

   double save_dt = 1e-3;

   cout << "[*] case: " << case_name << endl;
   // loading original case
   first_blood *fb = new first_blood(case_folder + case_name);
   fb->run();
   fb->save_results(save_dt);

   return 0;
}
