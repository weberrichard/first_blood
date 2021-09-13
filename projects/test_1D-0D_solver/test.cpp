#include "../../source/first_blood.h"

using namespace std;

int main(int argc, char* argv[])
{
	// string case_name = "P045H";
    string case_name = "Ferreira_heart";
   // string case_name = "lumped_test";
   // string case_name = "moc_lumped_test1";
   // string case_name = "moc_lumped_test2";
   // string case_name = "moc_lumped_test3";
   // string case_name = "moc_lumped_test4";
   // string case_name = "moc_lumped_test_heart";
   //string case_name = "moc_lumped_test_heart2";
   // string case_name = "moc_test";
   
	string case_folder = "../../models/";

   first_blood *fb = new first_blood(case_folder + case_name);
   cout << "const ok" << endl;

   fb->run();

   cout << "run ok" << endl;
   fb->save_results();
   cout << "save ok" << endl;

   cout << endl << endl;
   return 0;
}