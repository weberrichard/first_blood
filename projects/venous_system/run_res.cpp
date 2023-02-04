#include "../../source/first_blood.h"
#include <string>

using namespace std;

vector<double> vfr_bc(vector<double> t, double q_input_f);

int main(int argc, char* argv[])
{
   // basic stuff
	string case_folder = "models/";
   string case_name = "Cerebral_res";
   double R0, exp;
   if(argc == 4)
   {
      case_name = argv[1];
      R0 = stod(argv[2],0);
      exp = stod(argv[3],0);
   }
   else
   {
      cout << " Wrong number of input arguments: " << argc-1 << " right: 3" << endl;
      exit(-1);
   }

   cout << "[*] case: " << case_name << endl;

   // loading original case
   first_blood *fb = new first_blood(case_folder + case_name);

   // setting the resistance values
   double r0 = 0.000654648;
   vector<string> R0_id{"A57x","A58x","A64x","A65x","A71x","A74x","A72x","A75x","A76x","A78x","A80x","A82x","A84x","A86x","A88x","A90x","A91x","A93x","A92x","A94x","A100x","A102x"};
   vector<double> radius{0.000648,0.000648,0.00079425,0.00079425,0.000627025,0.000627025,0.000627025,0.000627025,0.0007315,0.0007315,0.000313525,0.000313525,0.000627025,0.000627025,0.000668775,0.000668775,0.0007525,0.0009,0.0007525,0.0007525,0.00058525,0.00058525};
   for(int i=0; i<R0_id.size(); i++)
   {
      int idx = fb->lum[0]->edge_id_to_index(R0_id[i]);
      if(idx>-1)
      {
         fb->lum[0]->edges[idx]->parameter = pow(r0/radius[i],exp)*R0;
      }
      else
      {
         cout << " R0_id " << R0_id[i] << " is not found." << endl;
      }
   }

   // running the simu
   fb->run();

   // saving results
   double save_dt = 1e-3;
   fb->save_results();

   return 0;
}

