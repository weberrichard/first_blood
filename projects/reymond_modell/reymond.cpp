#include "../../source/first_blood.h"
#include <string>

using namespace std;

int main(int argc, char* argv[])
{
	string case_folder = "../../models/";

   vector<string> case_names;
   case_names.push_back("Halasz_P045_resistance");
   case_names.push_back("Reymond_103");

   double save_dt = 1e-3;

   for(int i=0; i<case_names.size(); i++)
   {
      cout << "[*] " << case_names[i] << ":" << endl;
      first_blood *fb = new first_blood(case_folder + case_names[i]);
      fb->run();
      cout << " run OK" << endl;
      fb->save_results(save_dt);
      cout << " save OK" << endl;
   }

   cout << endl << endl;
   return 0;
}