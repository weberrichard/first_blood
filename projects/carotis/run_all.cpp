#include "../../source/first_blood.h"
#include <string>

using namespace std;

vector<double> vfr_bc(vector<double> t, double q_input_f);

int main(int argc, char* argv[])
{
   // basic stuff
	string case_folder = "models/";
   vector<string> case_names{"Carotis_2_Basic","Carotis_2_CE","Carotis_2_CI","Carotis_2_EI","Carotis_2_Mur","Carotis_2_WK3"};
   // vector<string> case_names{"Carotis_2_Mur"};

   double save_dt = 1e-3;

   // running flexible case
   for(int i=0; i<case_names.size(); i++)
   {
      cout << "[*] case: " << case_names[i] << " FLEX" << endl;

      first_blood *fb = new first_blood(case_folder + case_names[i]);
      fb->run();
      fb->save_results(save_dt,case_names[i]+"_flex");
   }

   // loading rigid case
   for(int i=0; i<case_names.size(); i++)
   {
      cout << "[*] case: " << case_names[i] << " RIGID" << endl;

      first_blood *fb = new first_blood(case_folder + case_names[i]);

      for(int j=0; j<fb->number_of_moc; j++)
      {
         for(int k=0; k<fb->moc[j]->number_of_edges; k++)
         {
            fb->moc[j]->edges[k]->elasticity_spring *= 100.;
            fb->moc[j]->edges[k]->elasticity_voigt *= 100.;
         }
      }
      fb->run();
      fb->save_results(save_dt,case_names[i]+"_rigid");
   }

   return 0;
}
