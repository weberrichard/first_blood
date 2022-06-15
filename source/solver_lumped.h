
/*===================================================================*\
								     solver_lumped
									 ---------------

	 first_blood
	 R. Weber
	 git: 
\*==================================================================*/

#ifndef SOLVER_LUMPED_H
#define SOLVER_LUMPED_H

#include	"file_io.h"

#include "/usr/include/eigen3/Eigen/Eigen"
#include <iostream>
#include <fstream>
#include <algorithm>

using namespace Eigen;
using namespace std;

class solver_lumped
{
public:
	solver_lumped(string file_name, string folder_path);
	~solver_lumped();

	// general functions
	// giving initial conditions
	void initialization();

	// loading the CSV file
	void load_model();

	// solving the equations
	vector<vector<double> > solve_one_step(double dt, vector<vector<double> > coefs);
	
	// clear/setting do_save_memory vars
	void clear_save_memory();
	void set_save_memory(vector<string> edge_list, vector<string> node_list);

	// saving output vars to file
	void save_results();
	void save_results(string folder_name);
	void save_results(string folder_name, vector<string> edge_list, vector<string> node_list);
	void save_results(double dt, string folder_name);
	void save_results(double dt, string folder_name, vector<string> edge_list, vector<string> node_list);

	// save model to file
	void save_model(string model_name, string folder_name);
	void save_initials(string model_name, string folder_name);
	void load_initials();

	// name of the model
	string name;

	// path of the model
	string input_folder_path;

	// outer nodes, i.e. outer boundaries to other models
	vector<string> boundary_main_node;
	vector<string> boundary_model_node;

	// index of moc, index of edge in moc, index of node in lumped
	vector<vector<int> > boundary_indices;

	// time of the simulation
	vector<double> time;

	// setting general constants
	void set_constants(double g, double rho, double nu, double mmHg, double p0);

	// setting non-SI parameters
	void set_non_SI_parameters();

	// parameters of elastance function
	double elastance_max = 2.5;
	double elastance_min = 0.06;
	double heart_rate = 75.6; // from Charlton2019

private:
	// general constants
	double gravity; // [m/s2]
	double density; // [kg/m3]
	double kinematic_viscosity; // [m2/s]
	double mmHg_to_Pa = 133.3616; // [Pa/mmHg] for converting inputs from mmHg to Pa
	double atmospheric_pressure; // Pa

	// Eigen vars for linear solver
	MatrixXd A;
	VectorXd b;

	class node
	{
	public:
		// variables and properties of NODEs
		string name;
		// type of node: node, ground
		string type;
		// in and outgoing edge indicies
		vector<int> edge_in, edge_out;
		// actual pressure for calculations
		double p; // mmHg
		// for virtual nodes, only used if connected to elastance
		double y; // mmHg
		// saving field variables
		bool do_save_memory = true;
		// containing pressure at nodes in time for storing
		vector<double> pressure; //Pa
		// initial condition for pressure
		double pressure_initial; // Pa
		double pres_ini_non_SI; // mmHg
		// whether the pressure is prescribed with a ground, true means p=0
		bool is_ground;
		// if the node is an outer boundary, ie connected to an other model
		bool is_master_node = false;
	};

	class edge
	{
	public:
		// name of the edge
		string name;
		// name of the nodes at the beginning and at the end
		string node_name_start, node_name_end;
		// index of the nodes at the beginning and at the end
		int node_index_start, node_index_end;
		// type of edge
		string type;
		// 0: resistance, 1: capacity, 2: elastance (time-varying capacity) 3: inductor, 4: voltage, 5: diode
		int type_code;
		// coefficient of the edge, e.g. R, C, L
		double parameter; // SI
		double par_non_SI; // non-SI
		// actual volume flow rate for calculations
		double vfr; // ml/s
		// saving field variables
		bool do_save_memory = true;
		// volume flow rate in time for storing
		vector<double> volume_flow_rate; // m3/s
		// initial condition for volume flow rate
		double volume_flow_rate_initial; // m3/s
		double vfr_ini_non_SI; // ml/s
	};

	// building the network, finding indicies
	void build_system();

	// general elastance function
	double elastance(double t);
	double elastance_derived(double t);

public:
	// containing nodes and edges in a vector
	vector<node*> nodes;
	vector<edge*> edges;
	
	int edge_id_to_index(string edge_id);
	int node_id_to_index(string node_id);

	// size of vectors
	int number_of_nodes, number_of_edges, number_of_master, number_of_elastance;
};

#endif // SOLVER_LUMPED_H