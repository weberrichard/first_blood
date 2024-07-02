
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
#include    "statistics.h"

#include "/usr/include/eigen3/Eigen/Eigen"
#include <iostream>
#include <fstream>
#include <algorithm>

using namespace Eigen;
using namespace std;


enum LumpedType { PerifCoronary0D, Perif0D, Heart0D };

class D0_transport;//declaration because solver_lumped needs it

class solver_lumped
{
public:
	solver_lumped(string file_name, string folder_path);
	~solver_lumped();

	// general functions
	// giving initial conditions
	void initialization(double hr);

	// loading the CSV file
	void load_model();

	// new nonlinear solver in fb, function fills regarding coefs in Jac and in f
	void set_newton_size();
	void coefficients_newton(double t_act);
	void initialization_newton(double t_act);
	void substitute_newton(double t_act);
	// variables
	MatrixXd Jac;
	VectorXd x, f;

	// control functions
	void update_parameters(double t_act);
	void myogenic_control(double t_act);

	// OLD solving the linear equations
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

	// index of moc, index of edge in moc, direction of edge, index of node in lumped
	vector<vector<int> > boundary_indices;

	// time of the simulation
	vector<double> time;
	double time_period;

	// setting general constants
	void set_constants(double g, double rho, double nu, double mmHg, double p0);

	// setting non-SI parameters
	void set_non_SI_parameters();

	// parameters of elastance function
	double elastance_max_nom = 2.5;
	double elastance_min_nom = 0.06;
	double heart_rate;

	// coronary modelling parameters from Reymond2009
	double alpha_coronary = 0.;
	double beta_coronary = 0.;

	// constants of myogenic control
	bool do_myogenic = false;
	double tao = 20.; // time constant
	double G = .9; // gain
	double sat1 = .55; // saturation 1
	double sat2 = 2.; // saturation 2
	time_average *q_ave, *p_ave, *C_ave, *R_fact, *x_myo_ts; // time period average values
	int q_idx=0, p_idx=0, C_idx; // index for average values, which element's average
	double x_myo; // acting signal
	double q_ref=0., p_ref=0.;


	// RBC transport in 0D
	D0_transport* RBClum;
	bool do_lum_RBC_transport = false;
	double fi_init_RBC_lum = 0.;

	// Haemoglobin saturation transport in 0D
	D0_transport* HBsatlum;
	bool do_lum_HB_sat_transport = false;
	double init_HB_sat_lum = 0.;

	// Plasma O2 concentration transport in 0D
	D0_transport* PlasmaO2lum;
	bool do_lum_PlasmaO2_transport = false;
	double init_PlasmaO2_lum = 0.;

	//O2 transport function
    void O2transport(double v, double dt, double dx, int n, double fiStartNodePlasma, double fiEndNodePlasma, double fiStartNodeHB, double fiEndNodeHB);
    double dCO2_plasma(double CO2_plasma_old, double HBsat_old, double C_RBC);
    double HBsat_equilibrium(double PO2);
    double turn_source(double t);
    void save_tissueO2(string folder_name, const vector<double>& st, const vector<double>& time);

    //tissue O2 concentration vector and scalar
    vector<double> tissueO2v;
    vector<double> tissueO2v_old;
    vector<double> tissueO2_save;
    double tissueO2s;
    double tissueO2s_old;

    vector<double> dPlasmaO2;
    
    //tissue O2 concentration initial condition
    double init_tissueO2 = 2.466237942122186e-3;
    //init function for tissue O2
    void init_lum_tissueO2();

    //constant parameters of O2 transport, literature: Bing dissertation and article
    double Dc = 1.6e-9; // [m2/s] O2 diffusion coefficient in plasma/blood
    double fi_c = 0.011303; // porosity (capillaries) [-]
    double fi_t = 0.988697; // porosity (tissue) =1-fi_c [-]
    double alpha_b = 3.11e-5; // [1/mmHg] solubility of oxygen in blood
    double alpha_t = 3.95e-5; // [1/mmHg] solubility of oxygen in tissue
    double hc = 1.0e-6; // [m] wall thickness of capillary walls
    double S_V_c = 4.74e5; // [1/m] surface to voulme ratio in capillaries
    double kc = 4.2e-14; // [m2/mmHg/s]
    double Mmax = 2.0e-4; // [1/s] ????
    //double C50 = 2.6e-5; // [m3/m3]
    double C50 = 2.6e-5; // [m3/m3]
    double taoO2 = 0.08;//s

    //parameters of the haemoglobin saturetion curve
    double L_HBsat = 1.251; // [-]
    double k_HBsat = 0.0676; // [1/mmHg]
    double b_HBsat = -0.274; // [-]
    double m_HBsat = 17.71; // [mmHg]

    // a parameter for converting the numner of O2 molecules/m3 to m3/m3
    double Z = 3.73e-17; // m3/1

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
		double RBC_fi0Dn; // RBC concenrtation
		double HBsat_0Dn; //haemodlobid saturation [1]
		double PlasmaO2_0Dn; //plasma O2 concentration [m3/m3]

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
		vector<double> parameter; // SI
		vector<double> par_non_SI; // non-SI
		double parameter_factor = 1.; // dimensionless factor, multiplies the parameter
		// actual volume flow rate for calculations
		double vfr; // ml/s
		// saving field variables
		bool do_save_memory = true;
		// volume flow rate in time for storing
		vector<double> volume_flow_rate; // m3/s
		// initial condition for volume flow rate
		double volume_flow_rate_initial; // m3/s
		double vfr_ini_non_SI; // ml/s

		bool is_open = false;//for diodes only, needed for transport
	};

	// building the network, finding indicies
	void build_system();

	// general elastance function
	double elastance(double t);
	double elastance(double t, vector<double> par);
	double elastance_derived(double t, vector<double> par);

public:
	// containing nodes and edges in a vector
	vector<node*> nodes;
	vector<edge*> edges;
	
	int edge_id_to_index(string edge_id);
	int node_id_to_index(string node_id);

	// size of vectors
	int number_of_nodes, number_of_edges, number_of_master, number_of_elastance, number_of_moc;

	//changing cross-section for the rtansport in 0D ("virtual 1D")
	double delta_V(int edge_index, int node_index);
};

//determines the number of divison points for virtual 1D
int NX(double L,double dx, int maxN);

void Virt1DforLum(vector<double> &fi_old, vector<double> &fi, double v, double dt, double dx, int n, double fiStartNode, double fiEndNode);

// Every lumped model gets one for eash type of transport. This handles the transport of substances in 0D
class D0_transport {//every 0D model gets one of this class
public:
    TransportType TType;
    bool do_save_results = false;
    //for simple peripherals wirh 4 RLC circuits
    vector<double> fi_arteriole, fi_capillary, fi_venulare, fi_vein;
    vector<double> fi_old_arteriole, fi_old_capillary, fi_old_venulare, fi_old_vein;
    double dx_arteriole,dx_capillary, dx_venulare, dx_vein;
    double L_arteriole, L_capillary, L_venulare, L_vein;
    double A_arteriole, A_capillary, A_venulare, A_vein;//from file, it is A_average-dA_average because of the changing cross-section
    int nx_arteriole, nx_capillary, nx_venulare, nx_vein;
    LumpedType LType;


    //for the heart model
    double fi_RA, fi_RV, fi_LA, fi_LV;//right atrium, right ventricle, left atrium, left ventricle, only nodes
    double fi_old_RA, fi_old_RV, fi_old_LA, fi_old_LV;
    //pulmonary circulation (together with heart model)
    vector<double> fi_pul_art, fi_pul_vein;//pulmonary circulation, 2 virtual 1D and three nodes
    vector<double> fi_old_pul_art, fi_old_pul_vein;
    double L_pul_art, L_pul_vein;
    double A_pul_art, A_pul_vein;
    double fi_lung;
    int nx_pul_art, nx_pul_vein;
    double dx_pul_art, dx_pul_vein;

    //for sawing concentration in time

    //Perif0D
    vector<double> fi_arteriole_start, fi_arteriole_end;
    vector<double> fi_capillary_start, fi_capillary_end;
    vector<double> fi_venulare_start, fi_venulare_end;
    vector<double> fi_vein_start, fi_vein_end;

    //heart, pulmanory circualtion
    vector<double> fi_RA_save, fi_RV_save, fi_LA_save, fi_LV_save;
    vector<double> fi_pul_art_start, fi_pul_art_end, fi_pul_vein_start, fi_pul_vein_end;



    double ml_to_m3 = 1.0e-6;

        

    D0_transport(LumpedType LType, vector<string> sv, TransportType TType, double concentration_init);

    void update_fi(double dt, double& masterFi, solver_lumped& lum_mod, double fi_vena_cava);

    void UpdatePerifLumNode(int LumNodeIndex, double fiLeft, double fiRight, solver_lumped& lum_mod);
    void update_lung_fi(double& fi_lung, double fiLeft, double fiRight, solver_lumped& lum_mod);
    void prescribe_lung_fi(solver_lumped& lum_mod, TransportType TType);

    void initialization();
    void save_variables();
    void save_results(string fn, const vector<double>& time, string model_name);
    void save_vector(string fname, const vector<double>& st, const vector<double>& en, const vector<double>& time);
    void save_vector(string folder_name, const vector<double>& vect, const vector<double>& time);
    void set_save_memory();
    vector<double> linear_dist(double avg, double dist, int len);

};

#endif // SOLVER_LUMPED_H