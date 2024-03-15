/*===================================================================*\
										 solver_moc
									 ---------------

	 first_blood
	 R. Weber
	 git: 
\*==================================================================*/

#ifndef SOLVER_MOC_H
#define SOLVER_MOC_H

#include "file_io.h"
#include "moc_edge.h"
#include "moc_node.h"
#include "statistics.h"
#include "solver_lumped.h"

#include "/usr/include/eigen3/Eigen/Eigen"

#include <sys/stat.h> // mkdir
#include <algorithm>

using namespace Eigen;

class TransportNodeCl; //declared

class solver_moc
{
public:
	solver_moc(string a_name, string a_folder);
	~solver_moc();

	// name of the model
	string name;
	string input_folder_path;

	// type of edge: vis
	string type;

	// division points in overall in the moc model
	int sum_division_points;

	// outer nodes, i.e. outer boundaries to other models
	vector<string> boundary_model_node;
	vector<string> boundary_main_node;

	// inner nodes and edges
	vector<moc_edge*> edges;
	vector<moc_node*> nodes;

	// number of elements
	int number_of_nodes, number_of_edges;

	// index of upstream boundary for interpolation
	vector<int> index_upstream;
	// upper boundary p-t
	vector<vector<double> > time_upstream;
	vector<vector<double> > value_upstream; // SI in code
	vector<int> type_upstream; // 0: pressure mmHg in file, 1: volume flow rate ml/s in file
	vector<int> node_upstream; // which
	// number of which period is the simulation
	vector<int> period;
	
	vector<double> olufsen_def_const{2.e6,-2253.,8.65e4}; // default constants for olufsen model

	// giving initial conditions
	void initialization(double pressure_initial,int mat_type);
	// initialazing Newton's method for 1D/0D boundaries
	void initialization_newton(VectorXd &x, int N, int moc_edge_index, int edge_end);
	// substituting the results back to field variables
	void substitute_newton(int moc_edge_index, int edge_end, double t_act, double p, double q);

	// filling up nodes edge_in, edge_out
	void build_system();

	// handling the boundaries: upstream and inner BC nodes, edge_inner is -1 if not inner iteration
	void boundaries(int edge_idx, double t_act);

	// calculate initial timesteps
	void timesteps();
	// searching the min timestep with idx
	double min_time(int &idx);

	// solving one time step, giving back actual time step
	double solve_one_step();

	// interpolate, save field vars
	void post_process();
	
	// performing the backward calculation on edge_id elemment
	void backward_solver(string node_id);

	// loading the model from csv file
	void load_model();
	// loading the time-pressure curve from CSV
	void load_time_series(string file_name);
	void convert_time_series();

	// setting basic constants
	void set_constants(double g, double rho, double nu, double mmHg, double p0, double nu_p, double cfl);

	// setting the full model to solve
	void full_tree();

	// clear/setting do_save_memory vars
	void clear_save_memory();
	void set_save_memory(vector<string> edge_list, vector<string> node_list);

	// saving output vars
	void save_results();
	void save_results(string folder_name);
	void save_results(string folder_name, vector<string> edge_list, vector<string> node_list);
	void save_results(double dt, string folder_name);
	void save_results(double dt, string folder_name, vector<string> edge_list, vector<string> node_list);

	// saving the model
	void save_model(string model_name, string folder_name);
	void save_pt_series(string model_name, string folder_name);

	// saving last values for further initial conditions
	void save_initials(string model_name, string folder_name);
	void load_initials();

	// from id to index nodes
	int node_id_to_index(string node_id);	
	// from id to index edges
	int edge_id_to_index(string edge_id);
	// finding all correspondning nodes to edges
	vector<int> edge_to_node(vector<int> edge_idx);

	//RBC transport for nodes
	TransportNodeCl* RBC_node_transport;

private:

	const double pi=3.14159265359;

	double gravity; // [m/s2]
	double density; // [kg/m3]
	double kinematic_viscosity; // [m2/s]
	double mmHg_to_Pa; // [Pa/mmHg] for converting inputs from mmHg to Pa
	double atmospheric_pressure; // Pa
	double pressure_initial; // Pa
	double beta; // exponent for wave velocity

	// heart p-t diagram
	vector<string> pt_file_name;
	
};

//handles transport for moc_nodes----------------------------------
//one class for RBC transport...
class TransportNodeCl {//for 1D 
public:
    TransportType TType;

    TransportNodeCl(TransportType TType);

    void update_fi(double& fiNode, moc_node* node, const vector<moc_edge*>& edges);
    void update_master_fi(double& fiNode, moc_node* node, const vector<moc_edge*>& edges, solver_lumped& lum_mod);
};




#endif // SOLVER_MOC_H


// founding the upward edge from a node, returning the edge index
//int backward_tree(string node_id);
// determining the downward tree from a node for the forward_solver, returning the indicies of edges
//void forward_tree(string node_id);
// for containing the forward tree edges
//vector<int> forward_edges;
// for containing the forward tree nodes
//vector<int> forward_nodes;
//vector<int> unique(vector<int> x);