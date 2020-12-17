/*===================================================================*\
                                 Edge
                            ---------------                            

  Abstract class, because it contains pure virtual functions i.e.
  it cannot be instantiated. Multiple classes are derived (e.g.
  Artery, ... ).

  This calss involves general variables like start_node_name,
  start_boundary_type etc. that every Edge requires.

    first_blood
    R. Weber
    git: 
\*==================================================================*/

#ifndef EDGE_H
#define EDGE_H

#include <string>
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>

using namespace std;

class edge
{
public:
   edge(string a_name);
   ~edge();

   string name; // name or ID of the edge
   string node_name_start, node_name_end; // name of the nodes at the beginning and at the end
   int node_index_start, node_index_end; // index of the nodes at the beginning and at the end

   // parameters: name, short name for shorter formulas
   // set_short_parameters function matches these
   double length; // m
   double nominal_diameter; // m
   double nominal_thickness; // m
   double viscosity;  // Ns/m2
   double elasticity_spring; // Ns/m2
   double elasticity_voigt; // Ns/m2
   double resistance_start; // 1/ms
   double resistance_end; // 1/ms
   double geodetic_height_start=0.; //m
   double geodetic_height_end=0.; //m
   int division_points; // pc.

   double gravity = 9.806; // m/s2
   double density = 1050.; // m/s2
   double kinematic_viscosity = 3.e-6; // m2/s
   double atmospheric_pressure = 1.e5; // Pa
   double beta = 2.0; // exponent for wave velocity, -

   // cointaining field variables in time
   vector<double> pressure_start; // Pa
   vector<double> pressure_end; // Pa
   vector<double> velocity_start;  // m/s
   vector<double> velocity_end; // m/s
   vector<double> wave_velocity_start; // m/s
   vector<double> wave_velocity_end; // m/s
   vector<double> total_deformation_start; // - 
   vector<double> total_deformation_end; // -
   vector<double> damper_deformation_start; // 
   vector<double> damper_deformation_end; // -
   vector<double> diameter_start; // m
   vector<double> diameter_end; // m
   vector<double> volume_flow_rate_start; // m3/s
   vector<double> volume_flow_rate_end; // m3/s
   vector<double> mass_flow_rate_start; // kg/s
   vector<double> mass_flow_rate_end; // kg/s

   // printing input parameters to console
   void print_input();
   // setting initial condition to field variables and setting short parameters
   void initialization();
   // setting upstream pressure p[0], only in the case of upstream_boundary
   void set_pressure_upstream(double p_in);
   // calculating the field variables at the new time step level
   void solve();
   // calculating the new timesteps for each division point
   double new_timestep();
   // interpolating in time and space to new timestep level to equidistant mesh
   void interpolate(double dt_real);
   // saving start and end field variables to vectors in time
   void save_field_variables();
   // updatin every field variable (a,d,epsz, ...) from v and p
   void update_field_variables(double dt);

   //---------------------
   // BOUNDARIES
   // [*] junctions, inner points
   // to calculate the nodal pressure
   vector<double> boundary_start_coefficients(double dt);
   vector<double> boundary_end_coefficients(double dt);
   // to determine the velocity and pressure of the edge
   void boundary_start_variables(double dt, double p, double v);
   void boundary_end_variables(double dt, double p, double v);

   // [*] upstream boundary node
   // if the node is upstream pressure BC, this calculates the velocity
   double upstream_boundary(double dt, double p_in);

   // [*] periferia points 
   // for periferia points solving p=p0
   double boundary_periferia(double dt, double p_out);

private:
   // variables for calculations of new field variables
   vector<double> dt; // time step for each point
   vector<double> x, xp, xq; // coordinates for field variables
   vector<double> p, pp, pq; // pressure, Pa
   vector<double> v, vp, vq; // velocity, m/s
   vector<double> a; // wave velocity, m/s
   vector<double> epsz; // relative deformation, -
   vector<double> epsz2; // relative deformation of voigt element, -
   vector<double> d; // diameter, m
   vector<double> A; // cross section area, m2
   vector<double> h; // geodetic height, m

   // short notations
   double l, dn, sn, eta2, E1, E2, Rs, Re, g, rho, nu, dx, hs, he, an, p0;
   int nx;
   void set_short_parameters(); // matching the longer and shorter parameters

   // calculating temp variables for characteristics
   double JL(int i);
   double JL(double dt, double p, double v, double a, double epsz, double epsz2, double d, double h, double xp, double x);
   double JR(int i);
   double JR(double dt, double p, double v, double a, double epsz, double epsz2, double d, double h, double xp, double x);

   // to get L and R positions
   double boundary_start_position(double dt);
   double boundary_end_position(double dt);
};

#endif // EDGE_H
