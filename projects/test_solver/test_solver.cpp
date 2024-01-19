#include "../../source/first_blood.h"
#include <string>

using namespace std;

int main(int argc, char* argv[])
{
   // basic stuff
	string case_folder = "../../models/";
   //string case_folder = "../../models/";
   double heart_rate = 75.6;  // if there is a heart model   
   double period_time = 60./heart_rate;

   // tested models
   vector<string> case_names;
   // case_names.push_back("v_perif_p0");
   // case_names.push_back("p0_perif_v");
   // case_names.push_back("q_perif_p0");
   // case_names.push_back("p0_perif_q");
   // case_names.push_back("p_perif_p0");
   // case_names.push_back("p0_perif_p");
   // case_names.push_back("p_Re_perif");
   // case_names.push_back("p_Rs_perif");
   // case_names.push_back("p1_p0");
   // case_names.push_back("p_node_perif_p0");
   // case_names.push_back("p_join_perif_p0");
   // case_names.push_back("p_0D");
   // case_names.push_back("p_0D_perif");
   // case_names.push_back("p_0D_p37");
   
   // case_names.push_back("Abel");
   case_names.push_back("Abel_ref3_myogenic");

   int st = 0; // 0 ="MacCormack", 1 = "MoC"
   string sn = "0D"; // "MacCormack", "MoC"
   string sf = "results"; // results_maccormack, "results_moc"

   cout << "[O] SOLVER: " << sn << endl;

   for(int i=0; i<case_names.size(); i++)
   {
      cout << " [*] case: " << case_names[i] << endl;
      // loading original case
      first_blood *fb = new first_blood(case_folder + case_names[i]);
      cout << "   + load: OK" << endl;

      //fb->time_end = sim_time;
      fb->material_type = 1; // setting to olufsen 1, linear 0
      //vector<double> olufsen_def_const{2.e6,-2253.,8.65e4}; // default constants for olufsen model
      //fb->material_const = olufsen_def_const;
      fb->time_period = period_time;
      fb->is_periodic_run = false;
      fb->solver_type = st;
      //fb->time_node = "n1";
    
      // running the simulation
      bool is_run_ok = fb->run();
      cout << "   + run: OK" << endl;

      if(is_run_ok)
      {
         //fb->save_results(sf + "/" + case_names[i]);
         fb->save_results();
         cout << "   + save: OK" << endl;
      }
   }

   return 0;
}
