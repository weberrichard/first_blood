/*===================================================================*\
										first_blood
									 ---------------

  Main first_blood class. Contains basic variables (e.g. vector for 
  node and edge) and functions (e.g. building the system topology for
  solving the equations).
 
	 first_blood
	 R. Weber
	 https://github.com/weberrichard/first_blood
\*==================================================================*/

#ifndef FIRST_BLOOD_H
#define FIRST_BLOOD_H

#include "file_io.h"
#include "solver_moc.h"
#include "solver_lumped.h"
#include "statistics.h"


#include "/usr/include/eigen3/Eigen/Eigen"

#include <string>
#include <iomanip>
#include <vector>
#include <cmath>
#include <iostream>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <stdio.h>

using namespace Eigen;

//transport stuff
//a first_blood object gets one of this class. This handles 1D transport for the moc edges
class Transport_1D {
public:
    TransportType TType;

    Transport_1D(TransportType TType);

    void update_fi(vector<double> v, vector<double>& fi, vector<double>& fi_new, double l, double dt, double fiStart, double fiEnd);
    void prescribe_node_fi_CO2(TransportType TType, double& finode);
};

class Transport_node {//for 1D nodes
public:
    TransportType TType;

    Transport_node(TransportType TType);

    void update_fi(double& fiNode, moc_node* node, const vector<moc_edge*>& edges);
    void update_master_fi(double& fiNode, moc_node* node, const vector<moc_edge*>& edges, solver_lumped& lum_mod, solver_moc& moc_mod);
};

class first_blood
{
public:
	first_blood(string folder_name);
	~first_blood();

	// vector of the models
	vector<solver_moc*> moc;
	vector<solver_lumped*> lum;
	vector<string> nodes;

	// run type of the simulation: forward, backward, lumped
	string run_type;

	// running the simulation for every moc and lumped
	// return true if succesful
	bool run();

	// solving lumped model as boundary condition
	void solve_lum(int index, double dt);
	void solve_lum_newton(int index, double dt);
	// finding lowest time step level
	double lowest_new_time(int &moc_idx, int &e_idx);

	// finding and setting the master nodes in every model
	void build_master();

	// initialization of fb object
	void initialization();

	// number of models
	int number_of_moc, number_of_lum, number_of_nodes;

	// file initialization from init folder
	bool init_from_file=false;

	// Constants for hydraulics
	double gravity = 9.806; // [m/s2]
	double density = 1055.; // [kg/m3]
	double kinematic_viscosity = 3.e-6; // [m2/s]
	double mmHg_to_Pa = 133.3616; // [Pa/mmHg] for converting inputs from mmHg to Pa
	double atmospheric_pressure = 1.e5; // Pa
	double pressure_initial = 0.*mmHg_to_Pa + atmospheric_pressure; // [Pa]
	double poisson_coefficient = 0.5;
	double courant_number = 0.9;
	int material_type=0; // 0: linear, 1: olufsen
	// vector<double> material_const; // actual used constants
	int solver_type = 0; // 0:maccormack, 1:moc with EE,

	// lumped time step if only lumped model exists
	double dt_lumped = 1.e-3;

	/// Loading the system from CSV
	bool load_ok;
	bool load_model();
	bool load_main_csv();

	/// Saving results to file
	void clear_save_memory(); // not saving anything to memory
	void set_save_memory(string model_name, string model_type, vector<string> edge_list, vector<string> node_list);
	double save_file_dt = 0.0; // time step of saving data in files, if 0, every data is saved
	void save_results(); // default folder name: case_name
	void save_results(string folder_name);
	void save_results(double dt); // default folder name: case_name
	void save_results(double dt, string folder_name);
	void save_results(string folder_name, string model_name, string model_type, vector<string> edge_list, vector<string> node_list); // saving specific time results to save time
	void save_results(double dt, string folder_name, string model_name, string model_type, vector<string> edge_list, vector<string> node_list); // saving specific time results to save time

	// Saving full model
	void save_model(string folder_name);
	void save_model(string model_name, string folder_name);

	// Saving last values for initial conditions for further simulations
	void save_initials(string model_name, string folder_name);
	void load_initials();

	// recording the timesteps
	double time_end;
	int time_counter;

	// for peridoic checking
	double time_end_min = 5.;
	double time_end_max = 500.;
	string time_node = "H";
	string time_var = "P";
	bool is_periodic_run = false;
	double time_val_old = -1.e10;
	bool is_run_end(double t_act, double t_old);

	// lum model id to index
	int lum_id_to_index(string lum_id);

	// path of the input file
	string input_folder_path;
	// name of the case without extension or folders
	string case_name;

	int period=0; // saving which period the calculation is
	//double heart_rate = 75.6; // from Charlton2019

	//for baroreflex
	 double L_B = 0.5936; //s
    double k_B = 0.1086; //1/mmHg
    double b_B = 0.6274; //s
    double x_B0 = 114.33; //mmHg
    //double x_B0 = 114.33; //mmHg
	double time_period;// = 60./heart_rate;
	string sys_moc;
	string sys_edge_name; // systole pressure of this edge
	moc_edge* sys_edge;
	bool do_baroreflex = false;
	double T_sum = 0.; //sum of completed cycles
	double T_last = 0.;
	void set_sys_edge_pointer();
	void init_time_periods_for_lum(double T);

	void O2_transport(int moc_idx, int si, int ei, int e_idx, double t_act);
	void CO2_transport(int moc_idx, int si, int ei, int e_idx, double t_act);
	void save_transport_var_for_lum(int moc_idx, int si, int ei);

	double baroreflex(double sys);
	int period_of_first_lum = 0;// which period the furthest lumped model is in


	// autoregulation stuff
	bool do_autoregulation = false;
	void autoregulation();
	
	// time averaged series
	time_average *map;


	void calculate_time_average();
	void save_time_average(string folder_name);
	void save_time_average(double dt, string folder_name);

	// RBC transport
	bool do_RBC_transport = false;
	Transport_1D* RBC1D; //class handling the 1D transport stuff
	TransportType TRBCType = RBC;
	Transport_node* RBC_node_transport;

	// Haemoglobin saturation
	bool do_HBsat_transport = false;
	Transport_1D* HBsat1D; //class handling the 1D transport stuff
	TransportType THBType = HB_O2_saturation;
	Transport_node* HB_O2_node_transport;

	// Plasma O2 concentration
	bool do_Plasma_O2_transport = false;
	Transport_1D* Plasma_O21D; //class handling the 1D transport stuff
	TransportType TPlasmaO2 = C_Plasma_O2;
	Transport_node* PlasmaO2_node_transport;

	//CO2 transport
	bool do_pla_CO2_transport = false;
    bool do_rbc_CO2_transport = false;
    bool do_pla_HCO3_transport = false;
    bool do_rbc_HCO3_transport = false;
    bool do_HbCO2_transport = false;

    Transport_node* transport_node_CO2_pla;
    Transport_node* transport_node_CO2_rbc;
    Transport_node* transport_node_HCO3_pla;
    Transport_node* transport_node_HCO3_rbc;
    Transport_node* transport_node_HbCO2;

    Transport_1D* CO2_pla_1D;
    Transport_1D* CO2_rbc_1D;
    Transport_1D* HCO3_pla_1D;
    Transport_1D* HCO3_rbc_1D;
    Transport_1D* HbCO2_1D;



	void check_valves();
	void connect_0D_edges();


private:
	// data of boundary for forward, and backward simulation
	class boundary
	{
	public:
		string node; // wich node
		string file_name; // filename where p-t or (v-t) series is
	};

public:
	boundary upstream_boundary;
};

#endif //FIRST_BLOOD_H 