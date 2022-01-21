#include "../../source/first_blood.h"
#include <string>

using namespace std;

vector<double> vfr_bc(vector<double> t, double q_input_f);

int main(int argc, char* argv[])
{
   // basic stuff
	string case_folder = "../../models/";
   string case_name = "Reymond_99_heart";
   //string case_name = "Reymond_103";
   double save_dt = 1e-3;
   double heart_rate = 75.6;  // if there is a heart model
   double period_time = 60./heart_rate;
   double sim_time = 50.*period_time;
   string heart_id = "heart_kim";

   cout << "[+] orig: ";

   // setting original case
   first_blood *fb = new first_blood(case_folder + case_name);
   fb->time_end = sim_time;

   // heart parameters
   int heart_index = fb->lum_id_to_index("heart_kim");
   fb->lum[heart_index]->heart_rate = heart_rate;
   fb->lum[heart_index]->edges[0]->parameter *= 1.0;
   fb->lum[heart_index]->edges[2]->parameter *= 1.0;
   fb->lum[heart_index]->edges[1]->parameter *= 2.0;
   fb->lum[heart_index]->edges[3]->parameter *= 2.0;

   // node resistance
   for(int i=0; i<fb->moc[0]->number_of_nodes; i++)
   {
      fb->moc[0]->nodes[i]->resistance = 25. * fb->moc[0]->nodes[i]->resistance;
   }

   // lum model parameters
   vector<string> perif_id_brain{"p25","p26","p27","p28","p29","p30","p31","p32","p33","p34","p35","p36","p37","p38","p39","p40","p41","p42","p43","p44","p45","p46"};
   vector<string> perif_id_hands{"p4","p5","p6","p15","p16","p17"};
   vector<string> perif_id_legs{"p7","p8","p9","p10","p11","p12","p13","p14"};
   vector<string> perif_id_spine{"p18","p19","p20","p21","p22","p23","p24"};
   // organizing some stuff for loops
   vector<vector<string> > perif_id_lum{perif_id_hands,perif_id_legs,perif_id_spine,perif_id_brain};
   vector<double> perif_res_orig{12.23,21.69,5.33,16.42};
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

   fb->clear_save_memory();
   string model_name = "arterial";
   string model_type = "moc";
   vector<string> el{"A44","A1","A12","A20","A5","A8","A46","A48"};
   vector<string> nl{"n8","n1","n6"};
   fb->set_save_memory(model_name,model_type,el,nl);

   fb->run();

   fb->save_results(save_dt,case_name,model_name,model_type,el,nl);

   // evaluate outputs 
   int n_out = 16;
   vector<double> output_orig(n_out);
   vector<double> PWV_length{0.61607,0.69615,0.21385,0.69524};

   // pressures
   int n8idx = fb->moc[0]->node_id_to_index("n8"); // radial artery
   output_orig[0] = (diastole(fb->moc[0]->nodes[n8idx]->pressure,fb->time,sim_time-period_time) - fb->atmospheric_pressure)/fb->mmHg_to_Pa;
   output_orig[1] = (systole(fb->moc[0]->nodes[n8idx]->pressure,fb->time,sim_time-period_time) - fb->atmospheric_pressure)/fb->mmHg_to_Pa;
   int n1idx = fb->moc[0]->node_id_to_index("n1"); // aorta
   output_orig[2] = (diastole(fb->moc[0]->nodes[n1idx]->pressure,fb->time,sim_time-period_time) - fb->atmospheric_pressure)/fb->mmHg_to_Pa;
   output_orig[3] = (systole(fb->moc[0]->nodes[n1idx]->pressure,fb->time,sim_time-period_time) - fb->atmospheric_pressure)/fb->mmHg_to_Pa;
   int n6idx = fb->moc[0]->node_id_to_index("n6"); // carotid
   output_orig[4] = (diastole(fb->moc[0]->nodes[n6idx]->pressure,fb->time,sim_time-period_time) - fb->atmospheric_pressure)/fb->mmHg_to_Pa;
   output_orig[5] = (systole(fb->moc[0]->nodes[n6idx]->pressure,fb->time,sim_time-period_time) - fb->atmospheric_pressure)/fb->mmHg_to_Pa;
   // volume flow rates
   int A44idx = fb->moc[0]->edge_id_to_index("A44"); // femoral artery
   output_orig[6] = average(fb->moc[0]->edges[A44idx]->volume_flow_rate_end,fb->time)*1.e6*60.;
   int A1idx = fb->moc[0]->edge_id_to_index("A1"); // cardiac output
   output_orig[7] = average(fb->moc[0]->edges[A1idx]->volume_flow_rate_start,fb->time)*1.e6*60.;
   int A12idx = fb->moc[0]->edge_id_to_index("A12"); // carotis
   output_orig[8] = average(fb->moc[0]->edges[A12idx]->volume_flow_rate_start,fb->time)*1.e6*60.;
   output_orig[9] = systole(fb->moc[0]->edges[A12idx]->volume_flow_rate_start,fb->time,sim_time-period_time)*1.e6*60.;
   int A20idx = fb->moc[0]->edge_id_to_index("A20"); // carotis
   output_orig[10] = average(fb->moc[0]->edges[A20idx]->volume_flow_rate_start,fb->time)*1.e6*60.;
   output_orig[11] = systole(fb->moc[0]->edges[A20idx]->volume_flow_rate_start,fb->time,sim_time-period_time)*1.e6*60.;

   // PWV
   output_orig[12] = PWV_length[0]/time_delay_correl(fb->moc[0]->edges[A1idx]->pressure_start,fb->moc[0]->edges[A44idx]->pressure_end, fb->time, sim_time-1.1*period_time);
   int A5idx = fb->moc[0]->edge_id_to_index("A5");
   output_orig[13] = PWV_length[1]/time_delay_correl(fb->moc[0]->edges[A5idx]->pressure_end,fb->moc[0]->edges[A44idx]->pressure_end, fb->time, sim_time-1.1*period_time);
   int A8idx = fb->moc[0]->edge_id_to_index("A8");
   //output_orig[14] = PWV_length[2]/time_delay_correl(fb->moc[0]->edges[A8idx]->velocity_start,fb->moc[0]->edges[A8idx]->velocity_end, fb->time, sim_time-1.1*period_time);
   output_orig[14] = (average(fb->moc[0]->edges[A8idx]->velocity_start,fb->time) + average(fb->moc[0]->edges[A8idx]->wave_velocity_start,fb->time) + average(fb->moc[0]->edges[A8idx]->velocity_end,fb->time) + average(fb->moc[0]->edges[A8idx]->wave_velocity_end,fb->time) )*.5;
   int A46idx = fb->moc[0]->edge_id_to_index("A46");
   int A48idx = fb->moc[0]->edge_id_to_index("A48");
   output_orig[15] = PWV_length[3]/time_delay_correl(fb->moc[0]->edges[A46idx]->pressure_start,fb->moc[0]->edges[A48idx]->pressure_end, fb->time, sim_time-1.1*period_time);

   cout << "OK" << endl;

   // starting sensitivity analysis
   // -----------------------------

   // setting original case
   first_blood *fb2 = new first_blood(case_folder + case_name);
   fb2->time_end = sim_time;

   int n_par = 48;
   double delta = 0.001;
   vector<vector<double> > par_dist(n_par,vector<double> (n_par,1.));
   vector<vector<double> > output_par(n_par,vector<double> (n_out,.0));

   //vector<double> pars(n_par,0.);
   for(int i=0; i<n_par; i++)
   {
      for(int j=0; j<n_par; j++)
      {
         if(i==j)
         {
            par_dist[i][j] = 1. + delta;
         }
      }
   }
   for(int i=0; i<n_par; i++)
   {
      cout << "[+] par " << i << "/" << n_par << ": ";
      // moc model edge IDs
      vector<string> moc_eid_brain{"A15","A17","A86","A85","A90","A89","A94","A93","A16","A20","A6","A5","A13","A84","A83","A12","A88","A87","A92","A91","A82","A75","A74","A73","A69","A78","A77","A76","A68","A70","A71","A72","A80","A56","A58","A57","A59","A65","A64","A61","A60","A81","A102","A67","A63","A103","A101","A100","A79","A66","A62"};
      vector<string> moc_eid_hands{"A8","A10","A11","A9","A7","A21","A23","A22","A24","A25"};
      vector<string> moc_eid_legs{"A42","A50","A51","A53","A52","A55","A54","A43","A45","A44","A47","A46","A48","A49"};
      vector<string> moc_eid_spine{"A4","A3","A19","A1","A95","A2","A14","A18","A26","A27","A29","A30","A31","A33","A32","A28","A34","A35","A36","A37","A38","A39","A40","A41"};
      vector<vector<string> > moc_eid{moc_eid_brain,moc_eid_hands,moc_eid_legs,moc_eid_spine};

      // moc model node IDs
      vector<string> moc_nid_brain{"n26","n27","n28","n29","n30","n32","n33","n34","n35","n36","n41","n42","n40","n39","n38","n37","n48","n43","n44","n45","n46","n47","n49","n50","n31"};
      vector<string> moc_nid_hands{"n8","n9","n11","n12"};
      vector<string> moc_nid_legs{"n14","n15","n16","n17","n18","n19"};
      vector<string> moc_nid_spine{"n1","n2","n3","n4","n6","n7","n10","n51","n52","n20","n21","n22","n23","n24","n25","n13"};
      vector<vector<string> > moc_nid{moc_nid_brain,moc_nid_hands,moc_nid_legs,moc_nid_spine};

      // heart parameters
      int heart_index2 = fb2->lum_id_to_index("heart_kim");
      fb2->lum[heart_index2]->heart_rate = heart_rate;
      for(int j=0; j<5; j++)
      {
         fb2->lum[heart_index2]->edges[j]->parameter = par_dist[i][j]*fb->lum[heart_index]->edges[j]->parameter;
         //pars[j] = fb->lum[heart_index]->edges[j]->parameter;
      }
      // elastance
      fb2->lum[heart_index2]->elastance_max = par_dist[i][5]*fb->lum[heart_index]->elastance_max;
      fb2->lum[heart_index2]->elastance_min = par_dist[i][6]*fb->lum[heart_index]->elastance_min;
      fb2->lum[heart_index2]->heart_rate = par_dist[i][7]*fb->lum[heart_index]->heart_rate;
      //pars[5] = fb->lum[heart_index]->edges[j]->elastance_max;
      //pars[6] = fb->lum[heart_index]->edges[j]->elastance_min;
      //pars[7] = fb->lum[heart_index]->edges[j]->heart_rate;

      // node resistance
      for(int j=0; j<moc_nid.size(); j++)
      {
         for(int k=0; k<moc_nid[j].size(); k++)
         {
            int idx = fb->moc[0]->node_id_to_index(moc_nid[j][k]);
            fb2->moc[0]->nodes[idx]->resistance = par_dist[i][8+j]*fb->moc[0]->nodes[idx]->resistance;
         }
      }

      // setting perif
      for(int j=0; j<perif_id_lum.size(); j++)
      {
         for(int k=0; k<perif_id_lum[j].size(); k++)
         {
            string id = perif_id_lum[j][k];
            int idx = fb->lum_id_to_index(id);
            fb2->lum[idx]->edges[0]->parameter = par_dist[i][12+j]*fb->lum[idx]->edges[0]->parameter;
            fb2->lum[idx]->edges[1]->parameter = par_dist[i][16+j]*fb->lum[idx]->edges[1]->parameter;
            fb2->lum[idx]->edges[2]->parameter = par_dist[i][20+j]*fb->lum[idx]->edges[2]->parameter;
         }
      }

      // edge parameters
      for(int j=0; j<moc_eid.size(); j++)
      {
         for(int k=0; k<moc_eid[j].size(); k++)
         {
            string id = moc_eid[j][k];
            int idx = fb->moc[0]->edge_id_to_index(id);
            fb2->moc[0]->edges[idx]->length = par_dist[i][24+j]*fb->moc[0]->edges[idx]->length; 
            fb2->moc[0]->edges[idx]->nominal_diameter_start = par_dist[i][28+j]*fb->moc[0]->edges[idx]->nominal_diameter_start; 
            fb2->moc[0]->edges[idx]->nominal_diameter_end = par_dist[i][28+j]*fb->moc[0]->edges[idx]->nominal_diameter_end; 
            fb2->moc[0]->edges[idx]->nominal_thickness_start = par_dist[i][32+j]*fb->moc[0]->edges[idx]->nominal_thickness_start; 
            fb2->moc[0]->edges[idx]->nominal_thickness_end = par_dist[i][32+j]*fb->moc[0]->edges[idx]->nominal_thickness_end; 
            fb2->moc[0]->edges[idx]->elasticity_spring = par_dist[i][36+j]*fb->moc[0]->edges[idx]->elasticity_spring;
            fb2->moc[0]->edges[idx]->elasticity_voigt = par_dist[i][40+j]*fb->moc[0]->edges[idx]->elasticity_voigt;
            fb2->moc[0]->edges[idx]->viscosity = par_dist[i][44+j]*fb->moc[0]->edges[idx]->viscosity;
         }
      }

      // setting vars to save
      fb2->clear_save_memory();
      fb2->set_save_memory(model_name,model_type,el,nl);

      // running the simulation
      fb2->run();

      // pressures
      output_par[i][0] = (diastole(fb2->moc[0]->nodes[n8idx]->pressure,fb2->time,sim_time-period_time) - fb2->atmospheric_pressure)/fb2->mmHg_to_Pa;
      output_par[i][1] = (systole(fb2->moc[0]->nodes[n8idx]->pressure,fb2->time,sim_time-period_time) - fb2->atmospheric_pressure)/fb2->mmHg_to_Pa;
      output_par[i][2] = (diastole(fb2->moc[0]->nodes[n1idx]->pressure,fb2->time,sim_time-period_time) - fb2->atmospheric_pressure)/fb2->mmHg_to_Pa;
      output_par[i][3] = (systole(fb2->moc[0]->nodes[n1idx]->pressure,fb2->time,sim_time-period_time) - fb2->atmospheric_pressure)/fb2->mmHg_to_Pa;
      output_par[i][4] = (diastole(fb2->moc[0]->nodes[n6idx]->pressure,fb2->time,sim_time-period_time) - fb2->atmospheric_pressure)/fb2->mmHg_to_Pa;
      output_par[i][5] = (systole(fb2->moc[0]->nodes[n6idx]->pressure,fb2->time,sim_time-period_time) - fb2->atmospheric_pressure)/fb2->mmHg_to_Pa;
      // volume flow rates
      output_par[i][6] = average(fb2->moc[0]->edges[A44idx]->volume_flow_rate_end,fb2->time)*1.e6*60.;
      output_par[i][7] = average(fb2->moc[0]->edges[A1idx]->volume_flow_rate_start,fb2->time)*1.e6*60.;
      output_par[i][8] = average(fb2->moc[0]->edges[A12idx]->volume_flow_rate_start,fb2->time)*1.e6*60.;
      output_par[i][9] = systole(fb2->moc[0]->edges[A12idx]->volume_flow_rate_start,fb2->time,sim_time-period_time)*1.e6*60.;
      output_par[i][10] = average(fb2->moc[0]->edges[A20idx]->volume_flow_rate_start,fb2->time)*1.e6*60.;
      output_par[i][11] = systole(fb2->moc[0]->edges[A20idx]->volume_flow_rate_start,fb2->time,sim_time-period_time)*1.e6*60.;
      // pulse wave velocity
      output_par[i][12] = PWV_length[0]/time_delay_correl(fb2->moc[0]->edges[A1idx]->pressure_start,fb2->moc[0]->edges[A44idx]->pressure_end, fb2->time, sim_time-1.1*period_time);

      output_par[i][13] = PWV_length[1]/time_delay_correl(fb2->moc[0]->edges[A5idx]->pressure_end,fb2->moc[0]->edges[A44idx]->pressure_end, fb2->time, sim_time-1.1*period_time);
      output_par[i][14] = (average(fb2->moc[0]->edges[A8idx]->velocity_start,fb2->time) + average(fb2->moc[0]->edges[A8idx]->wave_velocity_start,fb2->time) + average(fb2->moc[0]->edges[A8idx]->velocity_end,fb2->time) + average(fb2->moc[0]->edges[A8idx]->wave_velocity_end,fb2->time) )*.5;
      output_par[i][15] = PWV_length[3]/time_delay_correl(fb2->moc[0]->edges[A46idx]->pressure_start,fb2->moc[0]->edges[A48idx]->pressure_end, fb2->time, sim_time-1.1*period_time);

      cout << "OK" << endl;
   }

   string file_name = "output_par.csv";
   FILE *out_file;
   out_file = fopen(file_name.c_str(),"w");
   for(int i=0; i<output_par.size(); i++)
   {
      for(int j=0; j<output_par[i].size(); j++)
      {
         fprintf(out_file, "%9.7e,", output_par[i][j]);
      }
      fprintf(out_file, "\n");
   }
   fclose(out_file);

   file_name = "output_orig.csv";
   out_file = fopen(file_name.c_str(),"w");
   for(int i=0; i<output_orig.size(); i++)
   {
      fprintf(out_file, "%9.7e\n", output_orig[i]);
   }
   fclose(out_file);

   cout << endl;
   return 0;
}
