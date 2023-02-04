#include "../../source/first_blood.h"
#include <string>
#include <random>
#include <fstream>

using namespace std;

int main(int argc, char* argv[])
{

   // reference case and run case
   string case_folder = "../../models/";
   string case_ref_name = "Reymond_103";
   first_blood *fb_ref = new first_blood(case_folder + case_ref_name);
   double R1=0., R2=0., C=0.;
   for(int i=0; i<fb_ref->number_of_lum; i++)
   {
      R1 += fb_ref->lum[i]->edges[0]->parameter;
      R2 += fb_ref->lum[i]->edges[1]->parameter;
      C += fb_ref->lum[i]->edges[2]->parameter;
   }
   R1 /= (double)fb_ref->number_of_lum;
   R2 /= (double)fb_ref->number_of_lum;
   C /= (double)fb_ref->number_of_lum;

   double Rn=0.;
   for(int i=0; i<fb_ref->moc[0]->number_of_nodes; i++)
   {
      Rn += fb_ref->moc[0]->nodes[i]->resistance;
   }
   Rn /= (double)fb_ref->moc[0]->number_of_nodes;

   cout << "original: " << endl << " R1:" << R1 << "  R2: " << R2 << "  C: " << C << "   Rn: " << Rn << endl;

   // vpd parameters
   int n=4000, nok=0;
	case_folder = "vpd_scen2_3/";
   R1 = 0.; R2 = 0.; C = 0., Rn = 0.;
   for(int i=0; i<n; i++)
   {
      case_ref_name = "Reymond_103_2_vp_" + to_string(i);
      first_blood *fb = new first_blood(case_folder + case_ref_name);
      if(fb->load_ok == true)
      {
         nok++;
         for(int i=0; i<fb->number_of_lum; i++)
         {
            R1 += fb->lum[i]->edges[0]->parameter;
            R2 += fb->lum[i]->edges[1]->parameter;
            C += fb->lum[i]->edges[2]->parameter;
         }
         for(int i=0; i<fb->moc[0]->number_of_nodes; i++)
         {
            Rn += fb->moc[0]->nodes[i]->resistance;
         }
      }
   }
   R1 /= (double)fb_ref->number_of_lum;
   R2 /= (double)fb_ref->number_of_lum;
   C /= (double)fb_ref->number_of_lum;
   Rn /= (double)fb_ref->moc[0]->number_of_nodes;

   R1 /= (double)nok;
   R2 /= (double)nok;
   C /= (double)nok;
   Rn /= (double)nok;

   cout << "improved: " << endl << " R1:" << R1 << "  R2: " << R2 << "  C: " << C << "  Rn:" << Rn << endl;

   cout << endl << endl;
   return 0;
}
