#include "../../source/first_blood.h"
#include <string>

using namespace std;

int main(int argc, char* argv[])
{
	string case_folder = "../../models/";

   vector<string> case_names;
   case_names.push_back("artery_base");
   case_names.push_back("artery_resistance");
   case_names.push_back("artery_wind_kessel");

   for(int i=0; i<case_names.size(); i++)
   {
      cout << "[*] " << case_names[i] << ":";
      first_blood *fb = new first_blood(case_folder + case_names[i]);
      fb->run();
      fb->save_results();
      cout << " OK" << endl;
   }

   cout << endl << endl;
   return 0;
}