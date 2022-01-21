#include "../../source/first_blood.h"
#include <string>

using namespace std;

int main(int argc, char* argv[])
{
	string case_folder = "../../models/";

   string case_name = "Reymond_103";

   double save_dt = 1e-3;

   cout << "[*] " << case_name << ":" << endl;
   first_blood *fb = new first_blood(case_folder + case_name);
   fb->time_end = 100.;
   fb->clear_save_memory();
   string model_name = "arterial";
   string model_type = "moc";
   vector<string> el{"A1","A42","A50"};
   vector<string> nl;
   fb->set_save_memory(model_name,model_type,el,nl);
   fb->run();
   cout << " run OK" << endl;
   fb->save_results(case_name,model_name,model_type,el,nl);
   cout << " save OK" << endl;

   cout << endl << endl;
   return 0;
}