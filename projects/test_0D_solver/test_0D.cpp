#include "../../source/first_blood.h"

using namespace std;

int main(int argc, char* argv[])
{
	// string case_name = "P045H";
   string case_name = "Ferreira_heart";
   // string case_name = "lumped_test";
   
	string case_folder = "../../models/";

   first_blood *fb = new first_blood(case_folder + case_name);
   cout << "const ok" << endl;

   fb->run();

   for(int i=0; i<fb->number_of_lum; i++)
   {  
      for(int j=0; j<fb->lum[i]->number_of_nodes; j++)
      {
         cout << fb->lum[i]->nodes[j]->name  << endl;
      }
      for(int j=0; j<fb->lum[i]->number_of_edges; j++)
      {
         cout << fb->lum[i]->edges[j]->name << "  " << fb->lum[i]->edges[j]->parameter << endl;
      }
   }

   cout << "run ok" << endl;
   fb->save_results();
   cout << "save ok" << endl;

   cout << endl << endl;
   return 0;
}