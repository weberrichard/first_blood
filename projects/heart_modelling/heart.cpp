#include "../../source/first_blood.h"
#include <string>

using namespace std;

int main(int argc, char* argv[])
{
	string case_folder = "../../models/";
   string case_name = "Halasz_P045_heart";
   double save_dt = 1e-3;

   double hr = 50.8906;

   double R1,R2,L1,L2;
   if(argc == 5)
   {
      R1 = stod(argv[1],0); // 8e6
      L1 = stod(argv[2],0); // 6.7e4
      R2 = stod(argv[3],0); // 8e6
      L2 = stod(argv[4],0); // 6.9e4
   }
   else
   {
      cout << "Incorrect number of inputs: " << argc << endl << "Correct one is 5, exiting...." << endl;
      exit(-1);
   }

   // loading case
   first_blood *fb = new first_blood(case_folder + case_name);
   fb->lum[0]->heart_rate = hr;
   fb->time_end = 17.2*60./hr;

   // setting the values
   fb->lum[0]->edges[0]->parameter = R1;
   fb->lum[0]->edges[1]->parameter = L1;
   fb->lum[0]->edges[2]->parameter = R2;
   fb->lum[0]->edges[3]->parameter = L2;

   // actually running
   fb->run();

   string model_name = "heart_kim";
   string model_type = "lum";
   vector<string> nl{"aorta","left-ventricular"};
   vector<string> el;
   fb->save_results(save_dt, case_name, model_name, model_type, el, nl);

   return 0;
}