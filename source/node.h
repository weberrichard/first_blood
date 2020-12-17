/*===================================================================*\
                                  node
                            ---------------

  Node class for organizing every node property and function.
 
    first_blood
    R. Weber
    git: 
\*==================================================================*/

#ifndef NODE_H
#define NODE_H

#include <string>
#include <vector>
#include <cmath>

using namespace std;

class node
{
public:
   node(string a_name);
   ~node();

   // name of the node
   string name;

   // boundary type fo nodes
   // 0: junction (with or without resistance), 1: periferia, 2: heart
   string type;
   int type_code;

   // containing the INgoing and OUTgoing edges
   vector<int> edge_in;
   vector<int> edge_out;

   // pressure and volume_flow_rate in time
   vector<double> pressure; // in time, Pa
   vector<double> volume_flow_rate; // "leakage" volume flow rate, m3/s

   double resistance; // resistance coefficient, 1/ms
   double is_resistance; // if there is "leakage" this is 1, otherwise 0 and will eliminate this term
   double pressure_out; // outside pressure 
   double density = 1050.; // density of the fluid (blood)

   bool is_upstream_boundary = false; // for heart or upper boundary it is true

   // for intersection points;
   vector<double> boundary_coefficients();
   void boundary_variables(double p);

   // printing input parameters to console
   void print_input();

   // setting initial condition to field variables
   void initialization();

private:
   double R, Ri, p0, rho;
   void set_short_parameters(); // matching the longer and shorter parameters
   
};
#endif // NODE_H