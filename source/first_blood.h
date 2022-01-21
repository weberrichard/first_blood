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
#include "solver_moc.h"
#include "solver_lumped.h"
#include "statistics.h"

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

	// finding and setting the master nodes in every model
	void build_master();

	// initialization of fb object
	void initialization();

	// number of models
	int number_of_moc, number_of_lum, number_of_nodes;

	// Constants for hydraulics, note: there are constants in Edge.h
	double gravity = 9.806; // [m/s2]
	double density = 1055.; // [kg/m3]
	double kinematic_viscosity = 3.e-6; // [m2/s]
	double mmHg_to_Pa = 133.3616; // [Pa/mmHg] for converting inputs from mmHg to Pa
	double atmospheric_pressure = 1.e5; // Pa
	double pressure_initial = 100.*mmHg_to_Pa + atmospheric_pressure; // [Pa]
	double beta = 2.0; // exponent of wave velocity

	/// Loading the system from CSV
	bool load_ok;
	bool load_model();
	bool load_main_csv();

	/// Saving results to file
	void clear_save_memory(); // not saving anything to memory
	void set_save_memory(string model_name, string model_type, vector<string> edge_list, vector<string> node_list);
	double save_file_dt = 0.0; // time step of saving data in files, if 0, every data is saved
	void save_results(); // default folder name: case_name
	void save_results(double dt); // default folder name: case_name
	void save_results(string folder_name, string model_name, string model_type, vector<string> edge_list, vector<string> node_list); // saving specific time results to save time
	void save_results(double dt, string folder_name, string model_name, string model_type, vector<string> edge_list, vector<string> node_list); // saving specific time results to save time

	// Saving full model
	void save_model(string folder_name);
	void save_model(string model_name, string folder_name);

	// recording the timesteps
	vector<double> time;
	double time_end;

	// lum model id to index
	int lum_id_to_index(string lum_id);

	// path of the input file
	string input_folder_path;
	// name of the case without extension or folders
	string case_name;

	int period=0; // saving which period the calculation is

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