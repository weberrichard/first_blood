#include "../../source/first_blood.h"
#include <string>

using namespace std;

vector<double> vfr_bc(vector<double> t, double q_input_f);

int main(int argc, char* argv[])
{
   // basic stuff
	string case_folder = "";
   string case_name = "Reymond_99_heart";
   double save_dt = 1e-3;
   double heart_rate = 75.6;  // if there is a heart model
   double sim_time = 50.*60./heart_rate;
   string heart_id = "heart_kim";

   // handling inputs
   double res_brain_f, res_hands_f, res_legs_f, res_spine_f, res_heart_f;
   double cap_brain_f, cap_hands_f, cap_legs_f, cap_spine_f;
   double ind_heart_f;
   if(argc == 11)
   {
      res_brain_f = stod(argv[1],0);
      res_hands_f = stod(argv[2],0);
      res_legs_f = stod(argv[3],0);
      res_spine_f = stod(argv[4],0);
      cap_brain_f = stod(argv[5],0);
      cap_hands_f = stod(argv[6],0);
      cap_legs_f = stod(argv[7],0);
      cap_spine_f = stod(argv[8],0);
      res_heart_f = stod(argv[9],0);
      ind_heart_f = stod(argv[10],0);
   }
   else
   {
      cout << "Incorrect number of inputs (" << argc << "). Right one: 10" << endl;
      exit(-1);
   }

   // loading original case
   first_blood *fb = new first_blood(case_folder + case_name);
   fb->time_end = sim_time;
   int heart_index = fb->lum_id_to_index("heart_kim");
   fb->lum[heart_index]->heart_rate = heart_rate;
   // setting heart
   fb->lum[heart_index]->edges[0]->parameter *= res_heart_f;
   fb->lum[heart_index]->edges[2]->parameter *= res_heart_f;
   fb->lum[heart_index]->edges[1]->parameter *= ind_heart_f;
   fb->lum[heart_index]->edges[3]->parameter *= ind_heart_f;

   // setting the model parameters with argv values
   // moce node parameters
   for(int i=0; i<fb->moc[0]->number_of_nodes; i++)
   {
      fb->moc[0]->nodes[i]->resistance = 25. * fb->moc[0]->nodes[i]->resistance;
   }

   // lum model parameters
   vector<string> id_hands{"p4","p5","p6","p15","p16","p17"};
   vector<string> id_legs{"p7","p8","p9","p10","p11","p12","p13","p14"};
   vector<string> id_spine{"p18","p19","p20","p21","p22","p23","p24"};
   vector<string> id_brain{"p25","p26","p27","p28","p29","p30","p31","p32","p33","p34","p35","p36","p37","p38","p39","p40","p41","p42","p43","p44","p45","p46"};
   // organizing some stuff for loops
   vector<vector<string> > id_lum{id_hands,id_legs,id_spine,id_brain};
   vector<double> res_lum{res_hands_f,res_legs_f,res_spine_f,res_brain_f};
   vector<double> cap_lum{cap_hands_f,cap_legs_f,cap_spine_f,cap_brain_f};

   // setting perif
   for(int i=0; i<id_lum.size(); i++)
   {
      for(int j=0; j<id_lum[i].size(); j++)
      {
         string id = id_lum[i][j];
         int idx = fb->lum_id_to_index(id);
         fb->lum[idx]->edges[0]->parameter *= res_lum[i];
         fb->lum[idx]->edges[1]->parameter *= res_lum[i];
         fb->lum[idx]->edges[2]->parameter *= cap_lum[i];
      }
   }

   // fielad variable for saving to memory / files
   // variables from moc
   fb->clear_save_memory();
   string model_name = "arterial";
   string model_type = "moc";
   vector<string> el{"A44","A1","A12","A20","A5","A8","A46","A48"};
   vector<string> nl{"n8","n1","n6"};
   fb->set_save_memory(model_name,model_type,el,nl);

   // running the simulation
   fb->run();

   // saving stuff to file
   fb->save_results(save_dt,case_name,model_name,model_type,el,nl);
   //fb->save_results(save_dt);

   return 0;
}

//--------------------------------------------------------------
vector<double> vfr_bc(vector<double> t, double q_input_f)
{
   vector<double> q(t.size(),0.);

   // nominal values for q[ml/s] - t[s]
   // a0, a1, b1, a2, b2, a3, b3 ...
   // a0 + a1*cos(i*w*t) + b1*sin(i*w*t) + a2*cos(i*w*t) + ...
   vector<double> coef{104.2, 119.,142.7,-17.79,129.7,-46.54,49.81,-27.56,27.54,-32.46,16.4,-19.38,-2.409,-8.002,5.942,-15.87, 3.372};
   double w = 6.282;

   for(int i=0; i<coef.size(); i++)
   {  
      double c = coef[i]*q_input_f;
      
      if(i==0) // a0 term
      {
         for(int j=0; j<q.size(); j++)
         {
            q[j] += c;
         }
      }
      else if(i%2 == 1) // ai terms
      {
         for(int j=0; j<q.size(); j++)
         {
            q[j] += c*cos((i+1)/2*w*t[j]);
         }
      }
      else if(i%2 == 0) // bi terms
      {
         for(int j=0; j<q.size(); j++)
         {
            q[j] += c*sin(i/2*w*t[j]);
         }
      }
   }

   // converto to SI
   for(int i=0; i<q.size(); i++)
   {
      q[i] *= 1.e-6;
   }

   return q;
}