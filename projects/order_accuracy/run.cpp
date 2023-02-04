#include "../../source/first_blood.h"
#include <string>

using namespace std;

int main(int argc, char* argv[])
{
   // basic stuff
	string case_folder = "../../models/";
   string case_name = "Reymond_99_heart_ref3";
   double save_dt = 1e-3;
   double heart_rate = 75.6;  // if there is a heart model   
   double period_time = 60./heart_rate;
   double sim_time = 10.*period_time;
   bool is_periodic_run = true;
   bool init_from_file = false;

   clock_t begin_time;
   vector<double> p_aor_sys;
   vector<double> p_aor_dia;
   vector<double> car_out;
   vector<double> run_time;
   vector<int> nx_fact{1,2,4,8,16};

   int n_sim=nx_fact.size();
   for(int i=0; i<n_sim; i++)
   {
      cout << "case_idx: " << i << endl;
      // loading original case
      first_blood *fb = new first_blood(case_folder + case_name);
      //fb->time_end = sim_time;
      fb->time_period = period_time;
      fb->is_periodic_run = is_periodic_run;
      fb->init_from_file = init_from_file;

      // fielad variable for saving to memory / files
      // variables from moc
      fb->clear_save_memory();
      string model_name = "arterial";
      string model_type = "moc";
      vector<string> el{"A1","A5","A7","A8","A12","A52","A49"};
      vector<string> nl{};
      fb->set_save_memory(model_name,model_type,el,nl);

      // save from heart model
      string model_name2 = "heart_kim";
      string model_type2 = "lum";
      vector<string> el2{"D-mitral","D-aorta"};
      vector<string> nl2{"left-atrium","left-ventricular","aorta"};
      fb->set_save_memory(model_name2,model_type2,el2,nl2);

      for(int j=0; j<fb->moc[0]->number_of_edges; j++)
      {
         fb->moc[0]->edges[j]->division_points *= nx_fact[i];
      }

      begin_time = clock();
      // running the simulation
      bool is_run_ok = fb->run();
      double rt = double( clock () - begin_time ) /  CLOCKS_PER_SEC;
      rt /= fb->time_end;
      rt *= fb->time_period;
      run_time.push_back(rt);

      int A1idx = fb->moc[0]->edge_id_to_index("A1");
      double dia = diastole(fb->moc[0]->edges[A1idx]->pressure_start,fb->moc[0]->edges[A1idx]->time,fb->time_end-period_time);
      double sys = systole(fb->moc[0]->edges[A1idx]->pressure_end,fb->moc[0]->edges[A1idx]->time,fb->time_end-period_time);
      double ave = average(fb->moc[0]->edges[A1idx]->volume_flow_rate_start,fb->moc[0]->edges[A1idx]->time);

      p_aor_dia.push_back(dia);
      p_aor_sys.push_back(sys);
      car_out.push_back(ave);
      if(is_run_ok)
      {
         cout << "case_idx: " << i << " OK :)" << endl;
      }
      else
      {
         cout << "case_idx: " << i << " NOT ok :(" << endl;
      }

   }

   cout << "Results: " << endl;
   for(int i=0; i<n_sim; i++)
   {
      printf("%8.5e, %8.5e, %8.5e, %8.5e\n",p_aor_dia[i],p_aor_sys[i],car_out[i],run_time[i]);
   }
   cout << "Fertich! " << endl;

   return 0;
}


/*if(is_run_ok)
{
   // saving stuff to file
   //fb->save_results(save_dt);
   fb->save_results(save_dt,case_name,model_name,model_type,el,nl);
   fb->save_results(save_dt,case_name,model_name2,model_type2,el2,nl2);
   //fb->save_model("Reymond_99_heart_ref2","models");
   //fb->save_initials("Reymond_99_heart_ref","models");
}*/