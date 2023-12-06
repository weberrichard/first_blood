#include "../../source/first_blood.h"
#include <string>
#include <chrono>
#include <iostream>
#include <fstream>
using namespace std;

using std::chrono::high_resolution_clock;
using std::chrono::duration_cast;
using std::chrono::duration;
using std::chrono::milliseconds;

int main(int argc, char* argv[])
{
   // basic stuff
	//string case_folder = "../../models/test_cases/";
   string case_folder = "../../models/";
   double heart_rate = 75.6;  // if there is a heart model   
   double period_time = 60./heart_rate;
   ofstream fout;
   fout.open("runtimes_maccormack.txt", ios::app);
   //fout.open("runtimes_moc.txt");
   double save_dt = 1e-3;


   // tested models
   string case_name = argv[1];


   int st = 0; // 0, 1
   string sn = "MacCormack"; // "MacCormack", "MoC"
   string sf = "results_maccormack"; // results_maccormack, "results_moc"

   cout << "[O] SOLVER: " << sn << endl;

   for(int i=0; i<1; i++)
   {
      cout << " [*] case: " << case_name << endl;
      // loading original case
      first_blood *fb = new first_blood(case_folder + case_name);
      cout << "   + load: OK" << endl;

      // setting perif parameters
   vector<string> perif_brain{"p29","p30","p31","p32","p33","p34","p35","p36","p37","p38","p39","p40","p45","p46"};
   vector<string> perif_face{"p25","p26","p27","p28","p41","p42","p43","p44"};
   vector<string> perif_hands{"p4","p5","p6","p15","p16","p17"};
   vector<string> perif_legs{"p7","p8","p9","p10","p11","p12","p13","p14"};
   vector<string> perif_spine{"p18","p19","p20","p21","p22","p23","p24","p47"};
   // organizing some stuff for loops
   vector<vector<string> > id_lum{perif_brain, perif_face, perif_hands, perif_legs, perif_spine};
   // setting perif
   for(int i=0; i<id_lum.size(); i++)
   {
      for(int j=0; j<id_lum[i].size(); j++)
      {
         string id = id_lum[i][j];
         int idx = fb->lum_id_to_index(id);
      }
   }

   // setting artery parameters
   vector<string> edge_brain{"A5","A6","A20","A15","A12","A80","A72","A70","A71","A76","A78","A73","A74","A75","A82","A102","A16","A56","A57","A58","A59","A60","A64","A61","A65","A62","A63","A81","A79","A66","A67","A101","A103","A69","A68","A77","A100"};
   vector<string> edge_face{"A17","A85","A86","A89","A90","A93","A94","A13","A83","A84","A87","A88","A91","A92"};
   vector<string> edge_hands{"A3","A4","A7","A8","A9","A10","A11","A19","A21","A22","A23","A24","A25"};
   vector<string> edge_legs{"A42","A50","A51","A53","A52","A55","A54","A43","A44","A45","A47","A46","A48","A49"};
   vector<string> edge_spine{"A27","A29","A30","A31","A32","A33","A28","A34","A35","A36","A37","A38","A39","A40","A41","A1","A95","A2","A14","A18","A26"};
   vector<string> edge_coronary{"A96","A97","A98","A99"};
   // organizing some stuff for loops
   vector<vector<string> > id_edge{edge_brain, edge_face, edge_hands, edge_legs, edge_spine, edge_coronary};
   for(int i=0; i<id_edge.size(); i++)
   {
      for(int j=0; j<id_edge[i].size(); j++)
      {
         string id = id_edge[i][j];
         int idx = fb->moc[0]->edge_id_to_index(id);
      }
   }

   // setting the model parameters with argv values
   vector<string> node_brain{"n30","n31","n36","n37","n38","n39","n40","n41","n42","n43","n44","n45","n46","n47","n48","n49","n50"};
   vector<string> node_face{"n26","n27","n28","n29","n32","n33","n34","n35"};
   vector<string> node_hands{"n6","n7","n8","n9","n10","n11","n12"};
   vector<string> node_legs{"n14","n15","n16","n17","n18","n19"};
   vector<string> node_spine{"n20","n21","n22","n23","n24","n25","n51","n52","n13"};
   vector<string> node_coronary{"n53"};
   vector<vector<string> > id_node{node_brain, node_face, node_hands, node_legs, node_spine, node_coronary};
   for(int i=0; i<id_node.size(); i++)
   {  
      for(int j=0; j<id_node[i].size(); j++)
      {
         string id = id_node[i][j];
      }
   }

   fb->clear_save_memory();
   string model_name = "arterial";
   string model_type = "moc";
   vector<string> el{"A1","A2","A3","A5","A8","A12","A14","A18","A20","A27","A28","A35","A37","A39","A41","A43","A44","A46","A48","A95"};
   vector<string> nl{"H","n8","n1","n32"};
   fb->set_save_memory(model_name,model_type,el,nl);

      fb->time_period = period_time;
      fb->is_periodic_run = false;
      fb->solver_type = st;
      //fb->time_node = "n1";

      fb->material_type = 1; // olufsen model
      //vector<double> olufsen_def_const{2.e6,-2253.,8.65e4}; // default constants for olufsen model
      //fb->material_const = olufsen_def_const;
    
      // running the simulation
      auto t1 = high_resolution_clock::now();
      bool is_run_ok = fb->run();
      cout << "   + run: OK" << endl;
      auto t2 = high_resolution_clock::now();
      auto ms_int = duration_cast<milliseconds>(t2 - t1);
      fout << std::to_string(ms_int.count()) << endl;
      fout.close();
      
      vector<double> fitness_function(16,0.);
   if(is_run_ok)
   {
      // saving stuff to file
      //fb->save_results(save_dt);
      fb->save_results(save_dt,case_name,model_name,model_type,el,nl);
      //fb->save_model("Reymond_99_heart_ref2","models");
      //fb->save_initials("Reymond_99_heart_ref","models");

      // calculating fitness function
      int n8idx = fb->moc[0]->node_id_to_index("n8");
      int n1idx = fb->moc[0]->node_id_to_index("n1");
      int n32idx = fb->moc[0]->node_id_to_index("n32");
      fitness_function[0] = (diastole(fb->moc[0]->nodes[n8idx]->pressure,fb->moc[0]->nodes[n8idx]->time,fb->time_end-period_time) - fb->atmospheric_pressure)/fb->mmHg_to_Pa;
      fitness_function[1] = (systole(fb->moc[0]->nodes[n8idx]->pressure,fb->moc[0]->nodes[n8idx]->time,fb->time_end-period_time) - fb->atmospheric_pressure)/fb->mmHg_to_Pa;
      fitness_function[2] = (diastole(fb->moc[0]->nodes[n1idx]->pressure,fb->moc[0]->nodes[n1idx]->time,fb->time_end-period_time) - fb->atmospheric_pressure)/fb->mmHg_to_Pa;
      fitness_function[3] = (systole(fb->moc[0]->nodes[n1idx]->pressure,fb->moc[0]->nodes[n1idx]->time,fb->time_end-period_time) - fb->atmospheric_pressure)/fb->mmHg_to_Pa;
      fitness_function[4] = (diastole(fb->moc[0]->nodes[n32idx]->pressure,fb->moc[0]->nodes[n32idx]->time,fb->time_end-period_time) - fb->atmospheric_pressure)/fb->mmHg_to_Pa;
      fitness_function[5] = (systole(fb->moc[0]->nodes[n32idx]->pressure,fb->moc[0]->nodes[n32idx]->time,fb->time_end-period_time) - fb->atmospheric_pressure)/fb->mmHg_to_Pa;

      // volume flow rates
      int A44idx = fb->moc[0]->edge_id_to_index("A44");
      int A1idx = fb->moc[0]->edge_id_to_index("A1");
      int A12idx = fb->moc[0]->edge_id_to_index("A12");
      int A20idx = fb->moc[0]->edge_id_to_index("A20");
      fitness_function[6] = average(fb->moc[0]->edges[A44idx]->volume_flow_rate_end,fb->moc[0]->edges[A44idx]->time)*1.e6*60.;
      fitness_function[7] = average(fb->moc[0]->edges[A1idx]->volume_flow_rate_start,fb->moc[0]->edges[A1idx]->time)*1.e6*60.;
      fitness_function[8] = average(fb->moc[0]->edges[A12idx]->volume_flow_rate_start,fb->moc[0]->edges[A12idx]->time)*1.e6*60.;
      fitness_function[9] = systole(fb->moc[0]->edges[A12idx]->volume_flow_rate_start,fb->moc[0]->edges[A12idx]->time,fb->time_end-period_time)*1.e6*60.;
      fitness_function[10] = average(fb->moc[0]->edges[A20idx]->volume_flow_rate_start,fb->moc[0]->edges[A20idx]->time)*1.e6*60.;
      fitness_function[11] = systole(fb->moc[0]->edges[A20idx]->volume_flow_rate_start,fb->moc[0]->edges[A20idx]->time,fb->time_end-period_time)*1.e6*60.;

      // pulse wave velocity
      vector<string> pwv_aortic{"A1","A95","A2","A14","A18","A27","A28","A35","A37","A39","A41","A43","A44"};
      vector<string> pwv_car_fem{"A5","A3","A2","A14","A18","A27","A28","A35","A37","A39","A41","A43","A44"};
      vector<string> pwv_bra_rad{"A8"};
      vector<string> pwv_fem_ank{"A46","A48"};
      vector<vector<string> > id_pwv{pwv_aortic,pwv_car_fem,pwv_bra_rad,pwv_fem_ank};
      for(int i=0; i<id_pwv.size(); i++)
      {
         double lsum=0.;
         for(int j=0; j<id_pwv[i].size(); j++)
         {
            string id = id_pwv[i][j];
            int idx = fb->moc[0]->edge_id_to_index(id);
            double vs = average(fb->moc[0]->edges[idx]->velocity_start,fb->moc[0]->edges[idx]->time);
            double ve = average(fb->moc[0]->edges[idx]->velocity_end,fb->moc[0]->edges[idx]->time);
            double as = average(fb->moc[0]->edges[idx]->wave_speed_start,fb->moc[0]->edges[idx]->time);
            double ae = average(fb->moc[0]->edges[idx]->wave_speed_end,fb->moc[0]->edges[idx]->time);
            double l = fb->moc[0]->edges[idx]->length;
            lsum += l;
            fitness_function[12+i] += (vs+as + ve+ae)*.5 * l;
         }
         fitness_function[12+i] /= lsum;
      }
   }

   string file_name = "output.csv";
   FILE *out_file;
   out_file = fopen(file_name.c_str(),"w");
   for(int i=0; i<fitness_function.size(); i++)
   {
      fprintf(out_file, "%9.7e\n", fitness_function[i]);
   }
   fclose(out_file);
   }


   

   return 0;
}
