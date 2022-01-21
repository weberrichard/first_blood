#include "../../source/first_blood.h"
#include <string>
#include <random>
#include <fstream>

using namespace std;

int ts_find(vector<double> x, double x0);
vector<double> ts_crop(vector<double> x, int idx);
double ts_average(vector<double> t, vector<double> y);
double ts_min(vector<double> t, vector<double> y);
double ts_max(vector<double> t, vector<double> y);
vector<double> vfr_bc(vector<double> t, double sigma);

int main(int argc, char* argv[])
{
	string case_folder = "../../models/";
   string case_ref_name = "Reymond_103_2";

   int n_vp = 4000; // number of vp

   // parameters of distribution
   double mu_norm = 0.; // expected value of noise
   double sigma_norm_1 = 0.1; // length, thickness, damping
   double sigma_norm_2 = 0.2; // diameters, elasticity, q-t input

   double sigma_log = 0.472; // resistance, compliance
   double mu_log = -sigma_log*sigma_log*.5; // for having 0 expected value

   // random
   mt19937 engine;
   random_device dev{};
   engine.seed(dev());
   normal_distribution<double> norm_dist_1(mu_norm,sigma_norm_1);
   normal_distribution<double> norm_dist_2(mu_norm,sigma_norm_2);
   lognormal_distribution<double> log_dist(mu_log,sigma_log);

   // reference case and run case
   first_blood *fb_ref = new first_blood(case_folder + case_ref_name);
   first_blood *fb = new first_blood(case_folder + case_ref_name);

   // changing boundary type to volume_flow_rate
   fb->moc[0]->type_upstream = 1;
   fb->moc[0]->value_upstream.clear();
   fb->moc[0]->time_upstream.clear();
   double dt = 1e-3;
   for(int i=0; i<1000; i++)
   {
      fb->moc[0]->time_upstream.push_back(i*dt);
   }

   // for output values
   vector<double> p_rad_d, p_rad_s; // radial artery, dias and sys, pressure
   vector<double> p_aor_d, p_aor_s; // ascending aorta, dias and sys, pressure
   vector<double> p_car_d, p_car_s; // common carotid artery, dias and sys, pressure
   vector<double> q_fem_a_l, q_fem_a_r; // femoral artery, average, flow rate
   vector<string> cases_ok; // collecting the names of ok runs

   // nodal pressures
   int i_rad = 8; // 8 or 11, A9 start or A23 start
   int i_aor = 1; // 1, A95 start
   int i_car = 6; // 6, A3 end

   // volume flow rates of edges
   int i_fem_l = 44; // node 15 or 16, A50 end or A52 end
   int i_fem_r = 43; // node 18 or 19, A44 end or A46 end

   for(int i=0; i<n_vp; i++)
   {
      // name of vp
      string case_name = case_ref_name+"_vp_"+to_string(i);
      fb->case_name = case_name;

      // noise variable
      double noise;

      // distributing boundary conditions
      fb->moc[0]->value_upstream = vfr_bc(fb->moc[0]->time_upstream, sigma_norm_2);

      // distributing parameters
      // arterial
      for(int j=0; j<fb->moc[0]->number_of_edges; j++)
      {
         // length
         noise = norm_dist_1(engine);
         fb->moc[0]->edges[j]->length = (1. + noise) * fb_ref->moc[0]->edges[j]->length;

         // diameter at start
         noise = norm_dist_2(engine);
         fb->moc[0]->edges[j]->nominal_diameter_start = (1. + noise) * fb_ref->moc[0]->edges[j]->nominal_diameter_start;

         // diameter at end
         noise = norm_dist_2(engine);
         fb->moc[0]->edges[j]->nominal_diameter_end = (1. + noise) * fb_ref->moc[0]->edges[j]->nominal_diameter_end;

         // thickness at start
         noise = norm_dist_1(engine);
         fb->moc[0]->edges[j]->nominal_thickness_start = (1. + noise) * fb_ref->moc[0]->edges[j]->nominal_thickness_start;

         // thickness at end
         noise = norm_dist_1(engine);
         fb->moc[0]->edges[j]->nominal_thickness_end = (1. + noise) * fb_ref->moc[0]->edges[j]->nominal_thickness_end;

         // damping
         noise = norm_dist_1(engine);
         fb->moc[0]->edges[j]->viscosity = (1. + noise) * fb_ref->moc[0]->edges[j]->viscosity;

         // elasticity spring
         noise = norm_dist_2(engine);
         fb->moc[0]->edges[j]->elasticity_spring = (1. + noise) * fb_ref->moc[0]->edges[j]->elasticity_spring;

         // elasticity damper
         noise = norm_dist_2(engine);
         fb->moc[0]->edges[j]->elasticity_voigt = (1. + noise) * fb_ref->moc[0]->edges[j]->elasticity_voigt;
      }
      // WK at perif
      for(int j=0; j<fb->number_of_lum; j++)
      {
         for(int k=0; k<fb->lum[j]->number_of_edges; k++)
         {
            noise = log_dist(engine);
            fb->lum[j]->edges[k]->parameter = (1. + noise) * fb_ref->lum[j]->edges[k]->parameter;
         }
      }
      // arterial node leakage
      for(int j=0; j<fb->moc[0]->number_of_nodes; j++)
      {
         noise = log_dist(engine);
         fb->moc[0]->nodes[j]->resistance = (1. + noise) * fb_ref->moc[0]->nodes[j]->resistance;
      }

      // running simulations
      cout << "[*] start running vp_" << i << endl;
      bool is_run_ok = fb->run();
      if(is_run_ok)
      {
         cout << "[*] END   running vp_" << i << " ok" << endl;
      }
      else
      {
         cout << "[*] END   running vp_" << i << " NOT ok" << endl;
      }

      // saving output for tests
      //vector<string> edge_list, node_list;
      //edge_list.push_back("A8");
      //fb->save_results(1e-3, case_name,"arterial","moc",edge_list,node_list);
      //fb->save_results(1e-3);

      // saving model to file
      if(is_run_ok)
      {
         fb->save_model(case_name, "models/");
      }
      else
      {
         fb->save_model(case_name+"_NO", "models/");
      }

      // collecting outputs
      if(is_run_ok)
      {
         // saving the case name
         cases_ok.push_back(case_name);

         // index of start of the last period
         int index = ts_find(fb->time,4.);

         // saving rad dias and sys
         p_rad_d.push_back(ts_min(ts_crop(fb->time,index),ts_crop(fb->moc[0]->nodes[i_rad]->pressure,index)));
         p_rad_s.push_back(ts_max(ts_crop(fb->time,index),ts_crop(fb->moc[0]->nodes[i_rad]->pressure,index)));

         // saving aor dias and sys
         p_aor_d.push_back(ts_min(ts_crop(fb->time,index),ts_crop(fb->moc[0]->nodes[i_aor]->pressure,index)));
         p_aor_s.push_back(ts_max(ts_crop(fb->time,index),ts_crop(fb->moc[0]->nodes[i_aor]->pressure,index)));

         // saving aor dias and sys
         p_car_d.push_back(ts_min(ts_crop(fb->time,index),ts_crop(fb->moc[0]->nodes[i_car]->pressure,index)));
         p_car_s.push_back(ts_max(ts_crop(fb->time,index),ts_crop(fb->moc[0]->nodes[i_car]->pressure,index)));

         // saving fem q average
         q_fem_a_l.push_back(ts_average(ts_crop(fb->time,index),ts_crop(fb->moc[0]->edges[i_fem_l]->volume_flow_rate_end,index)));
         q_fem_a_r.push_back(ts_average(ts_crop(fb->time,index),ts_crop(fb->moc[0]->edges[i_fem_r]->volume_flow_rate_end,index)));
      }
      
      // writing outputs to file
      FILE *out_file;
      out_file = fopen("output.txt","w");
      for(int i=0; i<p_rad_d.size(); i++)
      {
         fprintf(out_file,"%s,%8.5e,%8.5e,%8.5e,%8.5e,%8.5e,%8.5e,%8.5e,%8.5e\n",cases_ok[i].c_str(),p_rad_d[i],p_rad_s[i],p_aor_d[i],p_aor_s[i],p_car_d[i],p_car_s[i],q_fem_a_l[i],q_fem_a_r[i]);
      }
      fclose(out_file);
   }


   cout << endl << endl;
   return 0;
}

//--------------------------------------------------------------
int ts_find(vector<double> x, double x0)
{
   int index=0;
   int i=0;
   while(i<x.size() && index==0)
   {
      if(x[i]<x0 && x[i+1]>x0)
      {
         index = i;
      }
      i++;
   }

   return index;
}

//--------------------------------------------------------------
vector<double> ts_crop(vector<double> x, int idx)
{
   vector<double> v;
   int n=x.size();
   for(int i=0; i<idx; i++)
   {
      v.push_back(x[n-idx+i]);      
   }

   return v;
}

//--------------------------------------------------------------
double ts_average(vector<double> t, vector<double> y)
{
   // average of a time series using trapezoid rule
   double m=0.;
   for(int i=0; i<t.size()-1; i++)
   {
      m += (t[i+1]-t[i])*(y[i+1]+y[i])*.5;
   }
   m /= (t.back()-t[0]);

   return m;
}

//--------------------------------------------------------------
double ts_min(vector<double> t, vector<double> y)
{
   // min (or diastole) of time series
   double min = y[0];
   for(int i=1; i<y.size(); i++)
   {
      if(y[i]<min)
      {
         min = y[i];
      }
   }

   return min;
}

//--------------------------------------------------------------
double ts_max(vector<double> t, vector<double> y)
{
   // max (or systole) of time series
   double max = y[0];
   for(int i=1; i<y.size(); i++)
   {
      if(y[i]>max)
      {
         max = y[i];
      }
   }

   return max;
}

//--------------------------------------------------------------
vector<double> vfr_bc(vector<double> t, double sigma)
{
   vector<double> q(t.size(),0.);

   // nominal values for q[ml/s] - t[s]
   // a0, a1, b1, a2, b2, a3, b3 ...
   // a0 + a1*cos(i*w*t) + b1*sin(i*w*t) + a2*cos(i*w*t) + ...
   vector<double> coef{104.2, 119.,142.7,-17.79,129.7,-46.54,49.81,-27.56,27.54,-32.46,16.4,-19.38,-2.409,-8.002,5.942,-15.87, 3.372};
   double w = 6.282;

   // random generating stuff
   mt19937 engine;
   random_device dev{};
   engine.seed(dev());
   // normal dist stuff
   double mu = 1.0;
   normal_distribution<double> distribution(mu,sigma);

   for(int i=0; i<coef.size(); i++)
   {  
      double noise = distribution(engine);
      double c = coef[i]*noise;
      
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
