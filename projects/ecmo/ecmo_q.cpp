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
   double sim_time = 1.3*period_time;
   bool init_from_file = false;

   // handling inputs
   if(argc == 4)
   {
      case_name = argv[1];
   }
   else
   {
      cout << "Incorrect number of inputs (" << argc << "). Right one: 4" << endl;
      exit(-1);
   }

   // loading original case
   first_blood *fb = new first_blood(case_folder + case_name);
   fb->clear_save_memory();

   // setting the heart parameters
   double heart_factor = stod(argv[2],0);
   int heart_index = fb->lum_id_to_index("heart_kim_lit");
   fb->lum[heart_index]->edges[3]->parameter[0]  *= heart_factor; // E_rv_max, E_max

   // setting the ecmo revolution number
   double rev_factor = stod(argv[3],0);
   string model_name2 = "heart_kim_lit"; // p8 for femfem and femcar, heart_kim_lit for venven
   string model_type2 = "lum";
   vector<string> el2{"V1"};
   vector<string> nl2{};
   if(rev_factor!=-1.0)
   {
      int lum_index = fb->lum_id_to_index("heart_kim_lit"); // p8 for femfem and femcar, heart_kim_lit for venven
      fb->lum[lum_index]->edges[5]->parameter_factor = rev_factor;
      fb->set_save_memory(model_name2,model_type2,el2,nl2);
   }

   string model_name = "arterial";
   string model_type = "moc";
   vector<string> el{"A1"};
   vector<string> nl{};
   fb->set_save_memory(model_name,model_type,el,nl);

   string model_name3 = "heart_kim_lit"; 
   string model_type3 = "lum";
   vector<string> el3;
   vector<string> nl3{"p_LA1", "p_RV2"};
   fb->set_save_memory(model_name3,model_type3,el3,nl3);

   fb->time_end = sim_time;

   // running the simulation
   bool is_run_ok = fb->run();

   string save_name = case_name + "_" + to_string(heart_factor);

   double co_ave, out=0., q_ecmo_ave;
   if(is_run_ok)
   {
      fb->save_results(save_name,model_name,model_type,el,nl);
      fb->save_results(save_name,model_name2,model_type2,el2,nl2);
      fb->save_results(save_name,model_name3,model_type3,el3,nl3);

      vector<double> co = fb->moc[0]->edges[0]->volume_flow_rate_start;
      vector<double> t = fb->moc[0]->edges[0]->time;
      int i_crop = crop_after_T(co,t,t.back()-period_time);
      vector<double> t2(t.begin()+i_crop,t.end());
      vector<double> co2(co.begin()+i_crop,co.end());
      co_ave = average(co2,t2);
      out += co_ave;

      if(rev_factor!=-1.0)
      {
         int lum_index = fb->lum_id_to_index("heart_kim_lit"); // p8 for femfem and femcar, heart_kim_lit for venven
         vector<double> q_ecmo = fb->lum[lum_index]->edges[5]->volume_flow_rate;
         vector<double> t = fb->lum[lum_index]->time;
         int i_crop = crop_after_T(q_ecmo,t,t.back()-period_time);
         vector<double> t2(t.begin()+i_crop,t.end());
         vector<double> co2(q_ecmo.begin()+i_crop,q_ecmo.end());
         q_ecmo_ave = average(co2,t2);
         out += q_ecmo_ave;
      }
   }

   std::ofstream file("result_" + to_string(heart_factor) + ".txt");
   file << std::setprecision(10) << out << "\n";
   file.close();

   return 0;
}
