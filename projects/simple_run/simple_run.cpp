#include "../../source/first_blood.h"
#include <string>

using namespace std;

int main(int argc, char* argv[])
{
   // basic stuff
	string case_folder = "../../models/";
   string case_name;
   double save_dt = 1e-3;
   double heart_rate = 75.6;  // if there is a heart model   
   double period_time = 60./heart_rate;
   double sim_time = 10.*period_time;
   bool init_from_file = false;

   // handling inputs
   if(argc == 2)
   {
      case_name = argv[1];
   }
   else
   {
      cout << "Incorrect number of inputs (" << argc << "). Right one: 1" << endl;
      exit(-1);
   }

   // loading original case
   first_blood *fb = new first_blood(case_folder + case_name);
   fb->time_end = sim_time;
   fb->time_period = period_time;
   fb->is_periodic_run = false;
   fb->init_from_file = init_from_file;
   //int heart_index = fb->lum_id_to_index("rats");
   //fb->lum[heart_index]->heart_rate = heart_rate;
   fb->heart_rate = heart_rate;

   // fielad variable for saving to memory / files
   // variables from moc
   //fb->clear_save_memory();
   //string model_name = "arterial";
   //string model_type = "moc";
   //vector<string> el{"A1","A5","A7","A8","A12","A52","A49","A16"};
   //vector<string> nl{};
   //fb->set_save_memory(model_name,model_type,el,nl);

   // save from heart model
   //string model_name2 = "heart_kim";
   //string model_type2 = "lum";
   //vector<string> el2{"D-mitral","D-aorta"};
   //vector<string> nl2{"left-atrium","left-ventricular","aorta"};
   //fb->set_save_memory(model_name2,model_type2,el2,nl2);

   // running the simulation
   bool is_run_ok = fb->run();

   if(is_run_ok)
   {
      // saving stuff to file
      fb->save_results(save_dt);
      //fb->save_results(save_dt,case_name,model_name,model_type,el,nl);
      //fb->save_results(save_dt,case_name,model_name2,model_type2,el2,nl2);
      //fb->save_model("Reymond_99_heart_ref2","models");
      //fb->save_initials("Reymond_99_heart_ref","models");
   }

   return 0;
}
