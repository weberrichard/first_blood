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

   srand((unsigned int) time(0));
   clock_t ido = clock();

   for(int I=0; I<case_name.size(); I++)
   {  
      cout << "[*] case: " << case_name[I] << endl;
      // loading original case
      first_blood *fb = new first_blood(case_folder + case_name[I]);
      fb->time_end = sim_time;
      int heart_index = fb->lum_id_to_index("heart_kim");
      fb->lum[heart_index]->heart_rate = heart_rate;
      fb->lum[heart_index]->edges[0]->parameter *= 1.1;
      fb->lum[heart_index]->edges[2]->parameter *= 1.1;

      // lum model parameters
      vector<string> perif_id_brain{"p25","p26","p27","p28","p29","p30","p31","p32","p33","p34","p35","p36","p37","p38","p39","p40","p41","p42","p43","p44","p45","p46"};
      vector<string> perif_id_hands{"p4","p5","p6","p15","p16","p17"};
      vector<string> perif_id_legs{"p7","p8","p9","p10","p11","p12","p13","p14"};
      vector<string> perif_id_spine{"p18","p19","p20","p21","p22","p23","p24"};
      // organizing some stuff for loops
      vector<vector<string> > perif_id_lum{perif_id_hands,perif_id_legs,perif_id_spine,perif_id_brain};
      vector<double> perif_res_orig{28.,30.,17.,27.};
      vector<double> perif_cap_orig{26.35,12.48,16.95,10.00};

      // setting perif
      for(int i=0; i<perif_id_lum.size(); i++)
      {
         for(int j=0; j<perif_id_lum[i].size(); j++)
         {
            string id = perif_id_lum[i][j];
            int idx = fb->lum_id_to_index(id);
            fb->lum[idx]->edges[0]->parameter *= perif_res_orig[i];
            fb->lum[idx]->edges[1]->parameter *= perif_res_orig[i];
            fb->lum[idx]->edges[2]->parameter *= perif_cap_orig[i];
         }
      }

      for(int i=0; i<fb->moc[0]->nodes.size(); i++)
      {
         fb->moc[0]->nodes[i]->resistance *= 35.;
      }

      // fielad variable for saving to memory / files
      // variables from moc
      //fb->clear_save_memory();
      //string model_name = "arterial";
      //string model_type = "moc";
      //vector<string> el{"A44","A1","A12","A20","A5","A8","A46","A48"};
      //vector<string> nl{"n8","n1","n6"};
      //fb->set_save_memory(model_name,model_type,el,nl);

      ido = clock();
      // running the simulation
      fb->run();
      cout << endl << "\n Run:   " << double(clock()-ido)/ CLOCKS_PER_SEC << " s" << endl;
      cout << " timesteps: " << fb->time.size() << endl;

      // saving stuff to file
      //fb->save_results(save_dt,case_name,model_name,model_type,el,nl);
      fb->save_results(save_dt);
   }

   return 0;
}
