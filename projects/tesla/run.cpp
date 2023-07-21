#include "../../source/first_blood.h"
#include <string>

using namespace std;

int main(int argc, char* argv[])
{
   // basic stuff
	string case_folder = "../../models/";
   double heart_rate = 75.6;  // if there is a heart model   
   double period_time = 60./heart_rate;

   // tested models
   vector<string> case_names;
   case_names.push_back("Tesla");
   case_names.push_back("Pipe");

   int st = 1; // 0, 1
   string sn = "MoC"; // "MacCormack", "MoC"
   string sf = "results"; // results_maccormack, "results_moc"

   cout << "[O] SOLVER: " << sn << endl;

   for(int i=0; i<case_names.size(); i++)
   {
      cout << " [*] case: " << case_names[i] << endl;
      // loading original case
      first_blood *fb = new first_blood(case_folder + case_names[i]);
      cout << "   + load: OK" << endl;

      //fb->time_end = sim_time;
      fb->material_type = 0; // setting to olufsen 1, linear 0
      // vector<double> olufsen_def_const{2.e6,-2253.,8.65e4}; // default constants for olufsen model
      // fb->material_const = olufsen_def_const;
      fb->time_period = period_time;
      fb->is_periodic_run = false;
      fb->solver_type = st;
      //fb->time_node = "n1";
    
      // running the simulation
      bool is_run_ok = fb->run();
      cout << "   + run: OK" << endl;

      if(is_run_ok)
      {
         fb->save_results();
         cout << "   + save: OK" << endl;
      }
   }

   return 0;
}
