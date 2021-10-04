#include "../../source/first_blood.h"
#include <string>

using namespace std;

int main(int argc, char* argv[])
{
	string case_folder = "models/";

   string case_name = argv[1];

   first_blood *fb = new first_blood(case_folder + case_name);
   fb->lum[0]->heart_rate = 70.588; // based on fft of ger_P035e
   fb->run();

   // saving specific results
   string save_folder = case_name;
   string model_to_save = "heart_kim";
   string model_type = "lum";
   vector<string> edge_list;
   edge_list.push_back("D-mitral");
   edge_list.push_back("D-aorta");
   vector<string> node_list;
   node_list.push_back("aorta");
   node_list.push_back("left-ventricular");

   fb->save_results(case_name, model_to_save, model_type, edge_list, node_list);

   return 0;
}