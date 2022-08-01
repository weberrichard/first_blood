#include "../../source/first_blood.h"
#include <string>

using namespace std;

vector<double> vfr_bc(vector<double> t, double q_input_f);

int main(int argc, char* argv[])
{
   // basic stuff
   double dt=1.e-3;
	string case_folder = "models/";

   string case_name;
   if(argc ==2)
   {
      case_name = argv[1];
   }
   else
   {
      cout << "One input argument is needed for case_name. Exiting..." << endl;
      exit(-1);
   }

   double heart_rate = 75.6;
   double period_time = 60./heart_rate;

   cout << "[*] case: " << case_name << endl;
   // loading original case
   first_blood *fb = new first_blood(case_folder + case_name);
   fb->time_period = period_time;
   //fb->is_periodic_run = true;

   // fielad variable for saving to memory / files
   // variables from moc
   fb->clear_save_memory();
   string model_name = "arterial";
   string model_type = "moc";
   vector<string> el{"A12","A16"};
   vector<string> nl{"H"};
   fb->set_save_memory(model_name,model_type,el,nl);

   fb->run();
   fb->save_results(dt);
   cout << "[*] case: " << case_name << " OK" << endl;

   return 0;
}
