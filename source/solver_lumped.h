
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
class D0_edge;

class solver_lumped
{
public:
	solver_lumped(string file_name, string folder_path);
	~solver_lumped();

	// general functions
	// giving initial conditions
	void initialization(double time_period);

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
	time_average  *p_ave; // time period average values
	int q_idx=0, p_idx=0, C_idx; // index for average values, which element's average
	double x_myo=0.; // acting signal
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
    void O2transport(double dt);
    void pulmonary_O2transport(double dt);
    double dCO2_plasma(double CO2_plasma_old, double HBsat_old, double C_RBC);
    double HBsat_equilibrium(double PO2);
    double turn_source(double t);
    void save_tissueO2(string folder_name, const vector<double>& time);
    vector<double> sin_2(double scale, int nx);
    void assign_perif_O2_params(vector<string> sv);
    void assign_haemogobin_sat_params(vector<string> sv);
    void assign_pulmonary_O2_params(vector<string> sv);

    //tissue O2 concentration vector and scalar
    vector<double> tissueO2v;
    vector<double> tissueO2_save;
    double tissueO2s;


    //CO2 transport
    void CO2transport(double dt);

    //tissue CO2 concentration vector and scalar
    vector<double> tissueCO2v;
    vector<double> tissueCO2_save;
    double tissueCO2s;

    
    //tissue O2 concentration initial condition
    double init_tissueO2 = 0.0266;
    //double init_tissueO2 = 2.2e-3;
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
    double Mmax = 2.5e-4; // [1/s] ????
    double C50 = 2.6e-5; // [m3/m3]
    double taoO2 = 0.5;//s

    //parameters of the haemoglobin saturetion curve
    double L_HBsat = 1.251; // [-]
    double k_HBsat = 0.0676; // [1/mmHg]
    double b_HBsat = -0.274; // [-]
    double m_HBsat = 17.71; // [mmHg]

    // a parameter for converting the numner of O2 molecules/m3 to m3/m3
    double Z = 3.73e-17; // m3/1

    //partial pressure of O2 in alveolars
    double PO2_alveolar = 100.; // [mmHg]
    double K_pul_O2 = 1.33e-7; // [m3/s/mmHg]
    double taoO2_p = 0.4;//s
    double K_pul_scale = 4.479e-4;

    //paramteres for metabolic response
    double x_met = 0.;
    double Ct_ref;
    double tao_met;
    double G_met;
    time_average *Ct_ave;
    bool do_metabolic_res = false;
    double sat1_met, sat2_met;
    void metabolic_response(double t_act);
    double vessel_dilation(int edgeindex);
    void set_0D_pointers();

    bool do_pul_O2_rtansport = false;
    bool do_per_O2_rtansport = false;


    //stored D0 capillary edges for oxygenation
    D0_edge* pul_cap_BH;
    D0_edge* pul_cap_RBC;
    D0_edge* pul_cap_PO2;

    D0_edge* per_cap_BH;
    D0_edge* per_cap_RBC;
    D0_edge* per_cap_PO2;

    //Baroreflex
    bool check_whole_period(double t_act);
    void update_period_time(double T_act_new);


    double T_act; //actual time period
    double T_last; //last time period
    double T_sum; //end of the last cycle
    int period = 0;


    //CO2 transport modelling
    //plasma CO2
    bool do_per_CO2_transport = false;
    
    D0_transport* CO2_pla_lum;
    bool do_lum_pla_CO2_transport = false;
    double fi_init_CO2_pla = 0.;

    //RBC cytoplasm CO2
    D0_transport* CO2_rbc_lum;
    bool do_lum_rbc_CO2_transport = false;
    double fi_init_CO2_rbc = 0.;

    //plasma HCO3
    D0_transport* HCO3_pla_lum;
    bool do_lum_pla_HCO3_transport = false;
    double fi_init_HCO3_pla = 0.;

    //RBC cytoplasm HCO3
    D0_transport* HCO3_rbc_lum;
    bool do_lum_rbc_HCO3_transport = false;
    double fi_init_HCO3_rbc = 0.;

    //CO2 with haemoglobin
    D0_transport* HbCO2_lum;
    bool do_lum_HbCO2_transport = false;
    double fi_init_HbCO2 = 0.;

    //capillary edges
 	D0_edge* per_cap_CO2_pla;
 	D0_edge* per_cap_CO2_rbc;
 	D0_edge* per_cap_HCO3_pla;
 	D0_edge* per_cap_HCO3_rbc;
 	D0_edge* per_cap_HbCO2;

 	//parameters
 	double tao_hco3_rbc_pla;
 	double tao_hco3_pla_rbc;
 	double tao_co2_pla_tis;
 	double tao_co2_tis_pla;
 	double tao_co2_rbc_pla;
 	double tao_co2_pla_rbc;

 	double tao_hco3;
 	double tao_hbco2;

 	double RQ;//respiratory quotient

 	double alpha_co2_rbc, alpha_co2_pla, alpha_hco3_rbc, alpha_hco3_pla, alpha_co2_tis;
 	double ksi_pla, ksi_rbc, ksi_c;
 	double fi_rbc, fi_pla;

 	//balance functions
 	double eta_hco3(double HCO3_rbc, double HCO3_pla, double CO2_pla, double CO2_rbc);
 	double eta_hb(double HbCO2, double CO2_pla, double CO2_rbc);

 	void assign_perif_CO2_params(vector<string> sv);
 	void init_lum_tissueCO2();

 	double init_tissueCO2 = 0.0266; //m3/m3


private:
	// general constants
	double gravity; // [m/s2]
	double density; // [kg/m3]
	double kinematic_viscosity; // [m2/s]
	double mmHg_to_Pa = 133.3616; // [Pa/mmHg] for converting inputs from mmHg to Pa
	double atmospheric_pressure; // Pa
	const double pi = 3.14159265359;
	const double ml_to_m3 = 1.0e-6;

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

		double CO2_pla_n; //CO2 concentration in plasma [m3/m3]
		double CO2_rbc_n; //CO2 concentration in RBC cytoplasm [m3/m3]
		double HCO3_pla_n; //HCO3 concentration in plasma [m3/m3]
		double HCO3_rbc_n; //HCO3 concentration in RBC cytoplasm [m3/m3]
		double HbCO2_n; //HbCO2 in concentration form [m3/m3]

		// whether the pressure is prescribed with a ground, true means p=0
		bool is_ground;
		// if the node is an outer boundary, ie connected to an other model
		bool is_master_node = false;

		//transport edges
		vector<D0_edge*> D0_edges_in_RBC = {};
		vector<D0_edge*> D0_edges_out_RBC = {};

		vector<D0_edge*> D0_edges_in_HBsat = {};
		vector<D0_edge*> D0_edges_out_HBsat = {};

		vector<D0_edge*> D0_edges_in_PlasmaO2 = {};
		vector<D0_edge*> D0_edges_out_PlasmaO2 = {};

		//for CO2 transport
		vector<D0_edge*> D0_edges_in_CO2_pla = {};
		vector<D0_edge*> D0_edges_out_CO2_pla = {};

		vector<D0_edge*> D0_edges_in_CO2_rbc = {};
		vector<D0_edge*> D0_edges_out_CO2_rbc = {};

		vector<D0_edge*> D0_edges_in_HCO3_pla = {};
		vector<D0_edge*> D0_edges_out_HCO3_pla = {};

		vector<D0_edge*> D0_edges_in_HCO3_rbc = {};
		vector<D0_edge*> D0_edges_out_HCO3_rbc = {};

		vector<D0_edge*> D0_edges_in_HbCO2 = {};
		vector<D0_edge*> D0_edges_out_HbCO2 = {};

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

		bool is_open = true;//for diodes only, needed for transport
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

	friend class D0_edge;
	void capillary_O2_transport(double dt);

	void autoregulation(double t_act);
	void update_R_fact();

	friend D0_transport;

};

//determines the number of divison points for virtual 1D
int NX(double L,double dx, int maxN);

void Virt1DforLum(vector<double> &fi, double v, double dt, double dx, int n, double fiStartNode, double fiEndNode);

//class for virtual 1D transport
class D0_edge{
public:
	bool is_diode = false;
	bool is_capacitor = false;
	bool is_elastance = false;


	double V0 = 0.;// heart chamber volume if p=0 for elastances or capacitances.

	D0_edge(string D0_name, double L, double A, int  nx, TransportType TType, double init, string node_s_name, string node_e_name, string corr_edge_name);
	string D0_name;
	double dx, L, A;
	int nx;
	bool is_per_capillary = false;//peripheral systemic capillary
	bool is_pul_capillary = false;//pulmonary capillary
	vector<double> fi;
	double init=0.;

	solver_lumped::node* node_start;
	solver_lumped::node* node_end;

	string node_s_name;
	string node_e_name;

	vector<double> fi_start;
	vector<double> fi_end;

	void update_edge();
	void save();

	TransportType TType;

	bool do_save_memory=false;

	//original edge, usually a resisto
	solver_lumped::edge* corr_edge;
	string corr_edge_name;


	void virt1D(double dt);
	const double ml_to_m3 = 1.0e-6;
	const double mmHg_to_Pa = 133.3616;

	void update_diode();
	void update_capacitor(double dt);
	void update_elastance(double dt, double E);
};


// Every lumped model gets one for eash type of transport. This handles the transport of substances in 0D
class D0_transport {//every 0D model gets one of this class
public:
    TransportType TType;
    bool do_save_results = false;
 
    double ml_to_m3 = 1.0e-6;
    bool do_tissue_transport = false;
    bool do_tissue_CO2_transport = false;



    D0_transport( TransportType TType);

    void update_fi(double dt, solver_lumped& lum_mod, double t_act);
    void prescribe_node_fi_CO2(TransportType TType, double& finode);

    //void UpdatePerifLumNode(int LumNodeIndex, double fiLeft, double fiRight, solver_lumped& lum_mod);
    void prescribe_node_fi(TransportType TType, double& finode);

    void save_variables();
    void save_results(string fn, const vector<double>& time, string model_name);
    void save_vector(string fname, const vector<double>& st, const vector<double>& en, const vector<double>& time);
    void save_vector(string folder_name, const vector<double>& vect, const vector<double>& time);
    void set_save_memory();
    vector<double> linear_dist(double avg, double dist, int len);


	void update_nodes(solver_lumped& lum_mod);
	void update_edges(double dt, solver_lumped& lum_mod, double t_act);
	void connect_0D_edges(solver_lumped& lum_mod);

    vector<D0_edge*> D0_edges;//virtual 1D elements
    const double mmHg_to_Pa = 133.3616;
};

#endif // SOLVER_LUMPED_H