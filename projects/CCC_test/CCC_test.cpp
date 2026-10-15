#include "../../source/first_blood.h"
#include <string>

using namespace std;

int main(int argc, char* argv[])
{
   // basic stuff
	string case_folder = "../../models/";
   double save_dt = 1e-2;
   //string case_folder = "../../models/";
   double heart_rate = 75.6;  // if there is a heart model

   double period_time = 60./heart_rate;

   // tested models
   vector<string> case_names;

   case_names.push_back("Abel_A62_63_control");


   int st = 0; // 0 ="MacCormack", 1 = "MoC"
   string sn = "MacCormack"; // "MacCormack", "MoC"
   string sf = "results"; // results_maccormack, "results_moc"

   cout << "[O] SOLVER: " << sn << endl;

   for(int i=0; i<case_names.size(); i++)
   {
      cout << " [*] case: " << case_names[i] << endl;
      // loading original case
      first_blood *fb = new first_blood(case_folder + case_names[i]);
      cout << "   + load: OK" << endl;
      //cout<< "alma";cout<< "alma";

      //fb->time_end = sim_time;
      fb->material_type = 1; // setting to olufsen 1, linear 0
      fb->is_periodic_run = false;
      fb->solver_type = st;


      fb->clear_save_memory();
      //
      string model_name1 = "arterial";
      string model_type1 = "moc";
      vector<string> el{"A1","A12","A16","A6","A20","A5","A15", "A70","A73"};
      vector<string> nl{};
      fb->set_save_memory(model_name1,model_type1,el,nl);

      
      vector<string> e3{};
      vector<string> n3{};
      
      vector<string> e20{"R1","R2","R3"};
      vector<string> n20{"Q1","P1"};
      //saving a lot of peripheral stuff for reference autoregualtion
      for (int k = 0; k < fb->lum.size(); k++) {
         for(int j = 4; j< 47; j++){
            if (fb->lum[k]->name == "p" + std::to_string(j)) {
              //fb->set_save_memory("p" + std::to_string(j),"CO2_transport++",e3,n3);
              fb->set_save_memory("p" + std::to_string(j),"lumped",e20,n20);
          }
       }
      }

    
      // running the simulation
      
      bool is_run_ok = fb->run();
      cout << "   + run: OK" << endl;

      if(is_run_ok)
      {
         //fb->save_results(sf + "/" + case_names[i]);
         //fb->save_results(save_dt);
         //cout << "   + save: OK" << endl;
         fb->save_results(case_names[i],model_name1 ,model_type1 ,el,nl);
        


         for (int k = 0; k < fb->lum.size(); k++) {
            for(int j = 4; j< 48; j++){
               if (fb->lum[k]->name == "p" + std::to_string(j)) {
                  //fb->save_results(case_names[i],"p" + std::to_string(j),"CO2_transport++",e3,n3);
                  fb->lum[k]->saveVectorsToFile(fb->lum[k]-> time, fb->lum[k]->nodes[5]->pressure, "results//" +case_names[i]+ "//p"+ std::to_string(j) +"//Q1.txt");
                  fb->lum[k]->saveVectorsToFile(fb->lum[k]-> time, fb->lum[k]->nodes[1]->pressure, "results//" +case_names[i]+ "//p"+ std::to_string(j) +"//P1.txt");
                  fb->lum[k]->saveVectorsToFile(fb->lum[k]-> time, fb->lum[k]->edges[1]->volume_flow_rate, "results//" +case_names[i]+ "//p"+ std::to_string(j) +"//R1.txt");
                  fb->lum[k]->saveVectorsToFile(fb->lum[k]-> time, fb->lum[k]->edges[2]->volume_flow_rate, "results//" +case_names[i]+ "//p"+ std::to_string(j) +"//R2.txt");
                  fb->lum[k]->saveVectorsToFile(fb->lum[k]-> time, fb->lum[k]->edges[3]->volume_flow_rate, "results//" +case_names[i]+ "//p"+ std::to_string(j) +"//R3.txt");
                  //fb->lum[k]->save_tissueO2("results//"+ case_names[i] +"//p" + std::to_string(j) + "//tissueO2.txt", fb->lum[k]-> time);
            }
            
         }
         }


      }

   }

   return 0;
}
