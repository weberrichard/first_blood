#include "../../source/first_blood.h"
#include <string>

using namespace std;

void save_output(first_blood& fb, const std::string& case_name, double save_dt, const std::vector<std::string>& sv);

int main(int argc, char* argv[])
{
   if (argc < 2) {
    cerr << "Usage: " << argv[0] << " <case_name>\n";
    return 1;
   }

	string case_folder = "../../models/";
   double save_dt = 0.;
 
   string case_name = argv[1];

   int st = 0; // 0 ="MacCormack", 1 = "MoC"
   string sn = "MacCormack"; // "MacCormack", "MoC"

   cout << "[O] SOLVER: " << sn << endl;

   cout << " [*] case: " << case_name << endl;
   first_blood fb(case_folder + case_name);
   cout << "   + load: OK" << endl;

   fb.material_type = 1; // setting to olufsen 1, linear 0

   fb.is_periodic_run = false;
   fb.solver_type = st;

   //setting what to save
   fb.clear_save_memory();

   std::ifstream file(case_folder + case_name + "/custom_output.txt");

   if (!file) {
      cout<< "Can't find" + case_folder + case_name + "/custom_output.txt"<<endl;
      return 1;
   }

   vector<vector<string>>save_requests{};

   std::string line;
   while (std::getline(file, line)) {
      line.erase(remove(line.begin(), line.end(), ' '), line.end());
      line.erase(remove(line.begin(), line.end(), '\n'), line.end());
      line.erase(remove(line.begin(), line.end(), '\r'), line.end());
      vector<string> sv = separate_line(line);

      //saving time, for interpolation
      if( (sv.size()==2) && (sv[0]=="save_dt") ){save_dt=stod(sv[1],0);}

      if(sv.size()<3){continue;}

      save_requests.push_back(sv); //saved for saving

      if( sv[1] == "moc" || sv[1] == "lum" ){
         vector<string> dummy{};
         vector<string> to_save(sv.begin() + 3, sv.end());

         if(sv[2] == "node"){
            fb.set_save_memory(sv[0], sv[1], dummy, to_save);
         }

         if(sv[2] == "edge"){
            fb.set_save_memory(sv[0], sv[1], to_save, dummy);
         }
      }

      if(sv[2] == "transport" && sv[1] == "lum" ){
         vector<string> dummy{};
         fb.set_save_memory(sv[0], "RBC_transport", dummy, dummy);
         fb.set_save_memory(sv[0], "HBsat_transport", dummy, dummy);
         fb.set_save_memory(sv[0], "PlasmaO2_transport", dummy, dummy);
         fb.set_save_memory(sv[0], "CO2_transport++", dummy, dummy);
      }
   }


   //run the simulation
   bool is_run_ok = fb.run();
   cout << "   + run: OK" << endl;

   //save outputs
   for (const auto& sv : save_requests) {
    save_output(fb, case_name, save_dt, sv);}

   return 0;
}


void save_output(first_blood& fb, const std::string& case_name, double save_dt, const std::vector<std::string>& sv)
{
    if (sv.size() < 3)
        return;

    std::vector<std::string> dummy;

    if (sv[1] == "moc" || sv[1] == "lum") {
        std::vector<std::string> to_save(sv.begin() + 3, sv.end());

        if (sv[2] == "node") {
            if (save_dt > 0.)
                fb.save_results(save_dt, case_name, sv[0], sv[1],
                                dummy, to_save);
            else
                fb.save_results(case_name, sv[0], sv[1],
                                dummy, to_save);
        }
        else if (sv[2] == "edge") {
            if (save_dt > 0.)
                fb.save_results(save_dt, case_name, sv[0], sv[1],
                                to_save, dummy);
            else
                fb.save_results(case_name, sv[0], sv[1],
                                to_save, dummy);
        }
    }

    if(sv[2] == "transport" && sv[1] == "lum" ){
         if (save_dt > 0.){
         fb.save_results(save_dt, case_name, sv[0], "RBC_transport", dummy, dummy);
         fb.save_results(save_dt, case_name, sv[0], "HBsat_transport", dummy, dummy);
         fb.save_results(save_dt, case_name, sv[0], "PlasmaO2_transport", dummy, dummy);
         fb.save_results(save_dt, case_name, sv[0], "CO2_transport++", dummy, dummy);}
         else{
         fb.save_results(case_name, sv[0], "RBC_transport", dummy, dummy);
         fb.save_results(case_name, sv[0], "HBsat_transport", dummy, dummy);
         fb.save_results(case_name, sv[0], "PlasmaO2_transport", dummy, dummy);
         fb.save_results(case_name, sv[0], "CO2_transport++", dummy, dummy);
         }
      }
}