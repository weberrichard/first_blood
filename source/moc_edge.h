/*===================================================================*\
										 moc_edge
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

#ifndef MOC_EDGE_H
#define MOC_EDGE_H

#include <string>
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>

#include "/usr/include/eigen3/Eigen/Eigen"

using namespace std;
using namespace Eigen;

class moc_edge
{
public:
	moc_edge(string a_ID);
	~moc_edge();

	string name; // name of the moc_edge
	string ID; // ID of the moc_edge
	string type; // type of the moc_edge: vis
	string node_name_start, node_name_end; // name of the nodes at the beginning and at the end
	int node_index_start, node_index_end; // index of the nodes at the beginning and at the end

	// parameters: name, short name for shorter formulas
	// set_short_parameters function matches these
	double length; // m
	double nominal_diameter_start; // m
	double nominal_diameter_end; // m
	double nominal_thickness_start; // m
	double nominal_thickness_end; // m
	double viscosity;  // Ns/m2
	double viscosity_factor;  // 0 or 1
	double elasticity_spring; // Ns/m2
	double elasticity_voigt; // Ns/m2
	double resistance_start; // 1/ms
	double resistance_end; // 1/ms
	double geodetic_height_start=0.; //m
	double geodetic_height_end=0.; //m
	int division_points; // pc.

	// for controlling the forward and backward calculations
	bool do_solve;

	double gravity; // m/s2
	double density; // m/s2
	double kinematic_viscosity; // m2/s
	double kinematic_viscosity_factor = 1.; // -
	double atmospheric_pressure; // Pa
	double beta; // exponent for wave velocity, -

	// saving field variables
	bool do_save_memory = true;

	// time step for inner iterations mainly
	double dt_act; // s
	vector<double> time; // s

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
	void print_vars();

	// setting initial condition to field variables and setting short parameters
	void initialization(double p_init);
	// setting upstream pressure p[0], only in the case of upstream_boundary
	void set_pressure_upstream(double p_in);

	//---------------------
	// FORWARD CALCULATION
	// calculating the field variables at the new time step level
	void solve();
	// calculating the new timesteps for each division point
	bool new_timestep();
	// interpolating in time and space to new timestep level to equidistant mesh
	void interpolate();
	void interpolate_hds();
	// saving start and end field variables to vectors in time
	void save_field_variables();
	// updatin every field variable (a,d,epsz, ...) from v and p at every point
	void update_variables();
	void save_initials(FILE* out_file);
	void set_initials(vector<double> ic);

	//---------------------
	// BOUNDARIES
	// master boundaries
	vector<double> boundary_master_start(double dt_master);
	vector<double> boundary_master_end(double dt_master);

	// new functions for Newton-Raphson
	vector<double> initialization_newton_start();
	vector<double> initialization_newton_end();
	vector<double> boundary_newton_start(double qp, double Ap, double pp, double t_act);
	vector<double> boundary_newton_end(double qp, double Ap, double pp, double t_act);

	// [*] junctions, inner points
	// newton iteration, start and end bc if necessary
	vector<MatrixXd> Jac;
	vector<VectorXd> y, f;
	void set_newton_size(int n1, int n2);

	// NEWEST junction handling with one equation NEWTON
	vector<double> junction_newton_start(double pp, double t_act);
	vector<double> junction_newton_end(double pp, double t_act);

	// OLD
	// to calculate the nodal pressure
	vector<double> boundary_start_coefficients(double t_act);
	vector<double> boundary_end_coefficients(double t_act);
	// to determine the velocity and pressure of the moc_edge
	void boundary_start_variables(double dt, double p, double v);
	void boundary_end_variables(double dt, double p, double v);

	// [*] upstream boundary node
	// if the node is upstream pressure BC, this calculates the velocity
	double upstream_boundary_p(double dt, double p_in);
	double downstream_boundary_p(double dt, double p_in);
	// if the node is upstream vfr BC, this calculates the pressure
	double upstream_boundary_q(double dt, double q_in);
	double downstream_boundary_q(double dt, double q_in);
	double upstream_boundary_v(double dt, double v_in, double &q);
	double downstream_boundary_v(double dt, double v_in, double &q);

	// [*] periferia points 
	// for periferia points solving p=p0
	double boundary_periferia_start(double dt, double p_out);
	double boundary_periferia_end(double dt, double p_out);

	//---------------------
	// BACKWARD CALCULATION
	// calculating backward in the moc_edge, returning [t,p_up]
	vector<vector<double> > backward_solver(vector<double> t_d, vector<double> p_d, vector<double> vfr_d);

//private:
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

	const double pi=3.14159265359;

	// for backward calculations
	double dt_back_max = .5e-3; // maximum timestep for backward calculation
	double dt_back; // real timestep for backward calculation
	int nt_back; // sizes of field variables
	vector<double> t_back; // equidistant time vector
	vector<double> tp_back; // location of new P points in time
	vector<double> x_back; // absolute location of x points
	vector<double> dx_back; // location of new P points in space with respect to previeous x

	void initialization_back(vector<double> time_in, vector<double> pressure_in, vector<double> velocity_in);
	// calculating x_min as the new space step for backward simulation
	double new_spacestep_back();
	// solving the characteristics at the end side from ith node
	void solve_back();
	// interpolating the end side of the moc_edge from ith node
	void interpolate_back(double dx_real);
	// reducing field var vectors
	void reduce_field_vectors();

	// short notations
	double l, dns, dne, sns, sne, eta2, eta2_f, E1, E2, Rs, Re, g, rho, nu, nu_f, dx, hs, he, ans, ane, p0, Ans, Ane;
	unsigned int nx; // number of division points
	void set_short_parameters(); // matching the longer and shorter parameters

	// calculating temp variables for characteristics
	double JL(int i);
	double JL(double dt, double p, double v, double a, double epsz, double epsz2, double d, double h, double xp, double x);
	double JR(int i);
	double JR(double dt, double p, double v, double a, double epsz, double epsz2, double d, double h, double xp, double x);
	// for backward calculation
	double JA(double dt, double p, double v, double a, double epsz, double epsz2, double d, double h, double xp, double x);
	double JB(double dt, double p, double v, double a, double epsz, double epsz2, double d, double h, double xp, double x);

	// to get L and R positions at boundaries
	double boundary_start_position(double dt);
	double boundary_end_position(double dt);

	// updating the ith field variables
	void update_ith_variables(int i, double ex, double p, double epsz2, double epsz);

	// for inner iteration for pp of periferia boundary
	double f_perif_start(double pp, double p_out, double dt, double &v_s);
	double f_perif_end(double pp, double p_out, double dt, double &v_e);
	double f_upstream_p(double pp, double p_in, double dt, double &v_s);
	double f_downstream_p(double pp, double p_in, double dt, double &v_e);
	double f_upstream_q(double pp, double q_in, double dt, double &v_s);
	double f_downstream_q(double pp, double q_in, double dt, double &v_e);
};

#endif // MOC_EDGE_H
