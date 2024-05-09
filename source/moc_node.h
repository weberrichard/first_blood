/*===================================================================*\
										 moc_node
									 ---------------

  moc_node class for organizing every moc_node property and function.
 
	 first_blood
	 R. Weber
	 git: 
\*==================================================================*/

#ifndef MOC_NODE_H
#define MOC_NODE_H

#include <string>
#include <vector>
#include <iostream>
#include <cmath>

using namespace std;


class moc_node
{
public:
	moc_node(string a_name);
	~moc_node();

	// name of the moc_node
	string name;

	// boundary type fo nodes
	// 0: junction (with or without resistance), 1: periferia, 2: heart
	string type;
	int type_code;

	// containing the INgoing and OUTgoing edges
	vector<int> edge_in;
	vector<int> edge_out;

	// saving field variables
	bool do_save_memory = true;

	// time vector
	vector<double> time;
	
	// pressure and volume_flow_rate in time
	vector<double> pressure; // in time, Pa
	vector<double> volume_flow_rate; // "leakage" volume flow rate, m3/s

	double resistance; // resistance coefficient, 1/ms
	double is_resistance; // if there is "leakage" this is 1, otherwise 0 and will eliminate this term
	double pressure_out; // outside pressure 
	double density; // density of the fluid (blood)
	double RBC_node_fi; // concenrtation of RBC

	// bool is_upstream_boundary = false;
	int upstream_boundary = -1; // for heart or upper boundary it is true
	bool is_master_node = false; // if it is connected to an other model, e.g lumped model
	int master_node_lum = -1; // index of lumped model if 

	// for intersection points;
	vector<double> boundary_coefficients();
	void boundary_variables(double p, double tact);

	// saving stuff to memory
	void save_field_variables(double t, double p, double q);

	// printing input parameters to console
	void print_input();

	// setting initial condition to field variables
	void initialization(double p_init, double RBC_init);

private:
	double R, Ri, p0, rho;
	void set_short_parameters(); // matching the longer and shorter parameters
	
};
#endif // MOC_NODE_H