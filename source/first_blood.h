/*===================================================================*\
                              first_blood
                            ---------------

  Main first_blood class. Contains basic variables (e.g. vector for 
  node and edge) and functions (e.g. building the system topology for
  solving the equations).
 
    staci3
    Cs. Hos, R. Weber, T. Huzsvar
    https://github.com/weberrichard/staci3
\*==================================================================*/

#ifndef FIRST_BLOOD_H
#define FIRST_BLOOD_H

#include "node.h"
#include "edge.h"

#include "/usr/include/Eigen/Eigen/Eigen"

#include <string>
#include <iomanip>
#include <vector>
#include <cmath>
#include <iostream>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <sys/stat.h> // mkdir

using namespace Eigen;

class first_blood
{
public:
   first_blood(string filename, double a_time_end);
   ~first_blood();

   // Basic Node and Edge list
   vector<node*> nodes;
   vector<edge*> edges;

   // Constants for hydraulics, note: there are constants in Edge.h
   double gravity = 9.81; // [m/s2]
   double density = 1050.; // [kg/m3]
   double kinematic_viscosity = 3e-6; // [m2/s]
   double mmHg_to_Pa = 133.3616; // [Pa/mmHg] for converting inputs from mmHg to Pa
   double atmospheric_pressure = 1.e5; // Pa
   
   // level of printing
   int printLevel = 0;

   /// Loading the system from CSV
   void load_system_csv();

   /// Saving results to file
   void save_results();

protected:
   // name of the case without extension or folders
   string case_name;
   // path of the input file
   string input_file_path;
   // path of the input folder
   string input_folder_path;

   int number_of_edges, number_of_nodes, number_of_timesteps;

   // recording the timesteps
   vector<double> time;

   // Creates the indicies for the nodes of edges and edges of nodes i.e. indexing the whole system and also filling up the edge vector
   void build_system();

   // Converting ID to index
   int node_id_to_index(string node_id);
   int edge_id_to_index(string edge_id);

   // containing the upstream boundary function, e.g. heart
   vector<double> time_upstream;
   vector<double> pressure_upstream;
   int index_upstream=0; // for registring the position of the interpolation
   int period=0; // saving which period the calculation is

   // length of the simulation
   double time_end;

private:
   // for helping the io_csv, it breaks a string to vector<string> separated by ,
   vector<string> separate_line(string line);

};

#endif //FIRST_BLOOD_H 