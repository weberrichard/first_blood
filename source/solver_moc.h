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

#include <sys/stat.h> // mkdir
#include <algorithm>

class solver_moc
{
public:
	solver_moc(string a_name, string a_folder);
	~solver_moc();

	// name of the model
	string name;
	string input_folder_path;

	// outer nodes, i.e. outer boundaries to other models
	vector<string> boundary_model_node;
	vector<string> boundary_main_node;

	// inner nodes and edges
	vector<moc_edge*> edges;
	vector<moc_node*> nodes;

	// number of elements
	int number_of_nodes, number_of_edges, number_of_timesteps;

	// index of upstream boundary for interpolation
	int index_upstream;
	// upper boundary p-t
	vector<double> time_upstream;
	vector<double> pressure_upstream;

	// number of which period is the simulation
	int period;

	// simulation time
	vector<double> time;

	// giving initial conditions
	void initialization(double pressure_initial);

	// filling up nodes edge_in, edge_out
	void build_system();

	// evaluating full calculation, several forward and backward calculation
	void full_solver(string node_id, double time_end);

	// perfroming the solver_moc forward
	// from a certain node_id to the perif
	void forward_solver(string node_id, double time_end);
	// handling the boundaries: upstream and inner BC nodes
	void boundaries(double dt);
	// interpolate, update fieldvars, save vars
	void post_process(double dt);

	// solving one time step, giving back actual time step
	double solve_one_step();

	// performing the backward calculation on edge_id elemment
	void backward_solver(string node_id);

	// loading the model from csv file
	void load_model();
	// loading the time-pressure curve from CSV
	void load_pt_series(string file_name);
	void convert_pt_series();

	// setting basic constants
	void set_constants(double g, double rho, double nu, double mmHg, double p0, double beta);

	// setting the full model to solve
	void full_tree();

	// saving output vars
	void save_results();
	void save_results(string folder_name);
	void save_results(string folder_name, vector<string> edge_list, vector<string> node_list);

	// from id to index nodes
	int node_id_to_index(string node_id);	
	// from id to index edges
	int edge_id_to_index(string edge_id);

private:
	// determining the downward tree from a node for the forward_solver, returning the indicies of edges
	void forward_tree(string node_id);
	// for containing the forward tree edges
	vector<int> forward_edges;
	// for containing the forward tree nodes
	vector<int> forward_nodes;
	vector<int> unique(vector<int> x);

	// founding the upward edge from a node, returning the edge index
	int backward_tree(string node_id);

	double gravity; // [m/s2]
	double density; // [kg/m3]
	double kinematic_viscosity; // [m2/s]
	double mmHg_to_Pa; // [Pa/mmHg] for converting inputs from mmHg to Pa
	double atmospheric_pressure; // Pa
	double pressure_initial; // Pa
	double beta; // exponent for wave velocity
	
};

#endif // SOLVER_MOC_H