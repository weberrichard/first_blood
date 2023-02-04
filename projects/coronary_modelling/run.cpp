#include "../../source/first_blood.h"
#include <string>

using namespace std;

vector<double> vfr_bc(vector<double> t, double q_input_f);

int main(int argc, char* argv[])
{
   // basic stuff
	string case_folder = "models/";
   string case_name = "Reymond_99_heart_ref3";
   double alpha,beta;
   if(argc == 4)
   {
      case_name = argv[1];
      alpha = stod(argv[2],0);
      beta = stod(argv[3],0);
   }
   else
   {
      cout << " Wrong number of argument: " << argc-1 << "   right: 3" << endl;
      cout << " Exiting..." << endl;
      exit(-1);
   }

   double save_dt = 1e-3;

   cout << "[*] case: " << case_name << endl;

   // loading original case
   first_blood *fb = new first_blood(case_folder + case_name);

   // setting coronary model
   vector<string> cor_perif{"p1","p2","p3"};
   for(int i=0; i<cor_perif.size(); i++)
   {
      int idx = fb->lum_id_to_index(cor_perif[i]);
      if(idx>-1)
      {
         fb->lum[idx]->alpha_coronary = alpha;
         fb->lum[idx]->beta_coronary = beta;
      }
   }

   // setting save
   fb->clear_save_memory();
   string model_name = "arterial";
   string model_type = "moc";
   vector<string> el{"A1","A96","A97","A98","A99"};
   vector<string> nl;
   fb->set_save_memory(model_name,model_type,el,nl);

   // actual run
   bool is_run_ok = fb->run();
   fb->save_results();

   // some feedback
   if(is_run_ok)
   {
      cout << "[*] case: " << case_name << " OK :)" << endl;
   }
   else
   {
      cout << "[*] case: " << case_name << " NOT OK :(" << endl;
   }

   return 0;
}
