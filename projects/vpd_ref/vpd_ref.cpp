#include "../../source/first_blood.h"
#include <string>

using namespace std;

vector<double> vfr_bc(vector<double> t, double q_input_f);

int main(int argc, char* argv[])
{
   // basic stuff
	string case_folder = "../../models/";
   string case_name = "Reymond_99_heart_2";
   double save_dt = 1e-3;
   double sim_time = 15.*1.0;
   double heart_rate = 60.;  // if there is a heart model

   // handling inputs
   double length_f, diameter_f, thickness_f;
   double res_brain_f, res_hands_f, res_legs_f, res_spine_f;
   double cap_brain_f, cap_hands_f, cap_legs_f, cap_spine_f;
   double damping_f, elasticity_1_f, elasticity_2_f;
   double res_nodes_f;
   double q_input_f;
   if(argc == 17)
   {
      length_f = stod(argv[1],0);
      diameter_f = stod(argv[2],0);
      thickness_f = stod(argv[3],0);
      res_brain_f = stod(argv[4],0);
      res_hands_f = stod(argv[5],0);
      res_legs_f = stod(argv[6],0);
      res_spine_f = stod(argv[7],0);
      cap_brain_f = stod(argv[8],0);
      cap_hands_f = stod(argv[9],0);
      cap_legs_f = stod(argv[10],0);
      cap_spine_f = stod(argv[11],0);
      damping_f = stod(argv[12],0);
      elasticity_1_f = stod(argv[13],0);
      elasticity_2_f = stod(argv[14],0);
      res_nodes_f = stod(argv[15],0);
      q_input_f = stod(argv[16],0);
   }
   else
   {
      cout << "Incorrect number of inputs (" << argc << "). Right one: 17" << endl;
      exit(-1);
   }

   // loading original case
   first_blood *fb = new first_blood(case_folder + case_name);
   fb->time_end = sim_time;
   int heart_index = fb->lum_id_to_index("heart_kim");
   fb->lum[heart_index]->heart_rate = heart_rate;

   // setting the model parameters with argv values
   // moc edge parameters
   for(int i=0; i<fb->moc[0]->number_of_edges; i++)
   {
      fb->moc[0]->edges[i]->length = length_f * fb->moc[0]->edges[i]->length;
      fb->moc[0]->edges[i]->nominal_diameter_start = diameter_f * fb->moc[0]->edges[i]->nominal_diameter_start;
      fb->moc[0]->edges[i]->nominal_diameter_end = diameter_f * fb->moc[0]->edges[i]->nominal_diameter_end;
      fb->moc[0]->edges[i]->nominal_thickness_start = thickness_f * fb->moc[0]->edges[i]->nominal_thickness_start;
      fb->moc[0]->edges[i]->nominal_thickness_end = thickness_f * fb->moc[0]->edges[i]->nominal_thickness_end;
      fb->moc[0]->edges[i]->viscosity = damping_f * fb->moc[0]->edges[i]->viscosity;
      fb->moc[0]->edges[i]->elasticity_spring = elasticity_1_f * fb->moc[0]->edges[i]->elasticity_spring;
      fb->moc[0]->edges[i]->elasticity_voigt = elasticity_2_f*fb->moc[0]->edges[i]->elasticity_voigt;
   }

   // moce node parameters
   for(int i=0; i<fb->moc[0]->number_of_nodes; i++)
   {
      fb->moc[0]->nodes[i]->resistance = res_nodes_f * fb->moc[0]->nodes[i]->resistance;
   }

   // changing boundary type to volume_flow_rate
   fb->moc[0]->type_upstream = 1;
   fb->moc[0]->value_upstream.clear();
   fb->moc[0]->time_upstream.clear();
   double dt = 1e-3;
   for(int i=0; i<1000; i++)
   {
      fb->moc[0]->time_upstream.push_back(i*dt);
   }
   fb->moc[0]->value_upstream = vfr_bc(fb->moc[0]->time_upstream, q_input_f);

   // lum model parameters
   vector<int> idx_hands{4,5,6,15,16,17};
   vector<int> idx_legs{7,8,9,10,11,12,13,14};
   vector<int> idx_spine{18,19,20,21,22,23,24};
   vector<int> idx_brain{25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46};
   // organizing some stuff for loops
   vector<vector<int> > idx_lum{idx_hands,idx_legs,idx_spine,idx_brain};
   vector<double> res_lum{res_hands_f,res_legs_f,res_spine_f,res_brain_f};
   vector<double> cap_lum{cap_hands_f,cap_legs_f,cap_spine_f,cap_brain_f};

   for(int i=0; i<idx_lum.size(); i++)
   {
      for(int j=0; j<idx_lum[i].size(); j++)
      {
         int idx = idx_lum[i][j];
         if(case_name == "Reymond_99" || case_name == "Reymond_99_heart_2")
         {
            idx -= 3;
         }
         fb->lum[idx]->edges[0]->parameter = res_lum[i] * fb->lum[idx]->edges[0]->parameter;
         fb->lum[idx]->edges[1]->parameter = res_lum[i] * fb->lum[idx]->edges[1]->parameter;
         fb->lum[idx]->edges[2]->parameter = cap_lum[i] * fb->lum[idx]->edges[2]->parameter;
      }
   }

   // fielad variable for saving to memory / files
   // variables from moc
   //fb->clear_save_memory();
   string model_name = "arterial";
   string model_type = "moc";
   vector<string> el{"A44","A1","A12","A20","A5","A8","A46","A48"};
   vector<string> nl{"n8","n1","n6"};
   fb->set_save_memory(model_name,model_type,el,nl);

   // running the simulation
   fb->run();

   // saving stuff to file
   //fb->save_results(save_dt,case_name,model_name,model_type,el,nl);
   fb->save_results(save_dt);

   // calculating distances for pwv: aortic, car-fem, bra-rad, fem-ank
   vector<double> distances(4,0.);
   // aortic pwv, aorta-femoral distance
   el = {"A1","A95","A2","A14","A18","A27","A28","A35","A37","A39","A41","A43","A44"};
   for(int i=0; i<el.size(); i++)
   {
      int idx = fb->moc[0]->edge_id_to_index(el[i]);
      distances[0] += fb->moc[0]->edges[idx]->length;
   }
   // carotis-femoral pwv
   el = {"A5","A3","A2","A14","A18","A27","A28","A35","A37","A39","A41","A43","A44"};
   for(int i=0; i<el.size(); i++)
   {
      int idx = fb->moc[0]->edge_id_to_index(el[i]);
      distances[1] += fb->moc[0]->edges[idx]->length;
   }
   // brachial-radial pwv
   el = {"A8"};
   for(int i=0; i<el.size(); i++)
   {
      int idx = fb->moc[0]->edge_id_to_index(el[i]);
      distances[2] += fb->moc[0]->edges[idx]->length;
   }
   // femoral-ankle pwv
   el = {"A46","A48"};
   for(int i=0; i<el.size(); i++)
   {
      int idx = fb->moc[0]->edge_id_to_index(el[i]);
      distances[3] += fb->moc[0]->edges[idx]->length;
   }

   ofstream out_file("distances.txt");
   for(int i=0; i<distances.size(); i++)
   {
      out_file << distances[i] << '\n';
   }
   out_file.close();

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