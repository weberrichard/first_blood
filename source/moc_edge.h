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
	double nominal_diameter_start,   nominal_diameter_end; // m
	double nominal_thickness_start,  nominal_thickness_end; // m
	double resistance_start,         resistance_end; // 1/ms
	double geodetic_height_start=0., geodetic_height_end=0.; //m
	double elasticity; // Ns/m2
	int division_points; // pc.

	// for controlling the forward and backward calculations
	bool do_solve;

	double gravity; // m/s2
	double density; // m/s2
	double kinematic_viscosity; // m2/s
	double kinematic_viscosity_factor = 1.; // -
	double atmospheric_pressure; // Pa
	double poisson_coefficient;
	int material_type=1; // 0: linear, 1: olufsen
	vector<double> material_const;
	double courant_number;

	// saving field variables
	bool do_save_memory = true;

	// time step for inner iterations mainly
	double dt_act; // s
	vector<double> time; // s

	// cointaining field variables in time
	vector<double> pressure_start,          pressure_end; // Pa
	vector<double> velocity_start,          velocity_end; // m/s
	vector<double> wave_speed_start,        wave_speed_end; // m/s
	vector<double> area_start,              area_end; // m2
	vector<double> volume_flow_rate_start,  volume_flow_rate_end; // m3/s
	vector<double> mass_flow_rate_start,    mass_flow_rate_end; // kg/s
	vector<double> RBC_concentration_start, RBC_concentration_end;// SI

	// printing input parameters to console
	void print_input();
	void print_vars();

	// setting initial condition to field variables and setting short parameters
	void initialization(double p_init, int mat_type);
	// setting upstream pressure p[0], only in the case of upstream_boundary
	void set_pressure_upstream(double p_in);

	//---------------------
	// FORWARD CALCULATION
	// calculating the field variables at the new time step level
	void solve_maccormack();
	void solve_moc();
	// calculating the new timesteps for each division point
	bool new_timestep();

	// riemann invariants
	double W1L(int j, double dt, double &J_L);
	double W2R(int j, double dt, double &J_R);

	// interpolated positions
	double right_position(int j, double dt);
	double left_position(int j, double dt);

	// saving start and end field variables to vectors in time
	void save_field_variables();
	// updatin every field variable (a,d,epsz, ...) from v and p at every point
	void update_variables();
	void save_initials(FILE* out_file);
	void set_initials(vector<double> ic);
	void update();

	// ---------------- \\
	//    BOUNDARIES    \\
	// ---------------- \\
	// to get L and R positions at boundaries
	double boundary_start_position(double dt);
	double boundary_end_position(double dt);

	// retunring the boundary variables to the edge class
	void boundary_substitute_start(double t_act, double p, double q);
	void boundary_substitute_end(double t_act, double p, double q);

	// new functions for Newton-Raphson
	double initialization_newton_start();
	double initialization_newton_end();
	vector<double> boundary_newton_start(double qp, double pp, double t_act);
	vector<double> boundary_newton_end(double qp, double pp, double t_act);

	// [*] junctions, inner points
	// newton iteration, start and end bc if necessary
	vector<MatrixXd> Jac;
	vector<VectorXd> y, f;
	void set_newton_size(int n1, int n2);

	// inner junction boundary condition
	vector<double> boundary_junction_start(double pp, double t_act);
	vector<double> boundary_junction_end(double pp, double t_act);

	// if the node is upstream pressure BC, this calculates the velocity
	double boundary_pressure_start(double dt, double p_in);
	double boundary_pressure_end(double dt, double p_in);
	double boundary_flowrate_start(double dt, double q_in);
	double boundary_flowrate_end(double dt, double q_in);
	double boundary_velocity_start(double dt, double v_in, double &q);
	double boundary_velocity_end(double dt, double v_in, double &q);

	// for inner iteration for pp of periferia boundary
	double f_pressure_start(double pp, double p_in, double dt, double &v_s);
	double f_pressure_end(double pp, double p_in, double dt, double &v_e);
	double f_flowrate_start(double pp, double q_in, double dt, double &v_s);
	double f_flowrate_end(double pp, double q_in, double dt, double &v_e);

	//get function(s) for transport
	vector<double> getVelocity();
	vector<double> getArea();

	vector<double> RBCfi, RBCfinew; //RBC concentration SI

private:
	// changing diameter along the vessel
	double nominal_area(double xp, double &dn_dx);
	double nominal_area(int i, double &dn_dx);
	double nominal_wall_thickness(double xp, double &sn_dx);
	double nominal_wall_thickness(int i, double &sn_dx);

	double wave_speed(double x, double A);
	double pressure(double x, double A, double &dp_dA, double &dp_dx, double &dF_dx);
	double area(double x, double p);

	// variables for calculations of new field variables
	vector<double> dt; // time step for each point
	vector<double> x; // coordinates for field variables
	vector<double> p, pnew; // pressure, Pa
	vector<double> v, vnew; // velocity, m/s
	vector<double> a, anew; // wave velocity, m/s
	vector<double> A, Anew; // diameter, m
	vector<double> h; // geodetic height, m

	const double pi=3.14159265359;

	// short notations
	double l, dns, dne, Ans, Ane, sns, sne, E, Rs, Re, g, rho, nu, nu_f, dx, hs, he, p0, nu_p, cfl, beta;
	double k1, k2, k3; // olufsen model const
	unsigned int nx; // number of division points
	void set_short_parameters(); // matching the longer and shorter parameters

	// calculating temp variables for characteristics
	double JL(int i);
	double JL(double dt, double p, double v, double a, double d, double xp);
	double JR(int i);
	double JR(double dt, double p, double v, double a, double d, double xp);

};

#endif // MOC_EDGE_H


// double f_perif_start(double pp, double p_out, double dt, double &v_s);
// double f_perif_end(double pp, double p_out, double dt, double &v_e);

// OLD
// to calculate the nodal pressure
// vector<double> boundary_start_coefficients(double t_act);
// vector<double> boundary_end_coefficients(double t_act);

// [*] periferia points 
// for periferia points solving p=p0
// double boundary_periferia_start(double dt, double p_out);
// double boundary_periferia_end(double dt, double p_out);

// updating the ith field variables
// void update_ith_variables(int i, double p);

// vector<double> boundary_master_start(double dt_master);
// vector<double> boundary_master_end(double dt_master);

// for backward calculations
/*double dt_back_max = .5e-3; // maximum timestep for backward calculation
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
// for backward calculation
double JA(double dt, double p, double v, double a, double d, double xp);
double JB(double dt, double p, double v, double a, double d, double xp);

// calculating backward in the moc_edge, returning [t,p_up]
vector<vector<double> > backward_solver(vector<double> t_d, vector<double> p_d, vector<double> vfr_d);
*/

// interpolating in time and space to new timestep level to equidistant mesh
//void interpolate();
//void interpolate_hds();