#include "moc_node.h"

using namespace std;

moc_node::moc_node(string a_name)
{
   name = a_name;
}

//--------------------------------------------------------------
moc_node::~moc_node(){}

//--------------------------------------------------------------
void moc_node::print_input()
{
	printf("\n %8s, %8s, %6.4f, %6.4f, %6.4f, %6.4f", type.c_str(), name.c_str(), 0., p0, R, Ri);
}

//--------------------------------------------------------------
void moc_node::initialization(double p_init)
{
	// clearing time variables
   pressure.clear();
   volume_flow_rate.clear();

   // setting short versions of parameters
   set_short_parameters();

   // saving initial conditions
   if(type_code == 0 || type_code == 2)
   {
   	pressure.push_back(p_init);
   }
   else
   {
   	pressure.push_back(pressure_out);
   }
   volume_flow_rate.push_back(0.);

   is_upstream_boundary = false;
   is_master_node = false;
}

//--------------------------------------------------------------
void moc_node::set_short_parameters()
{
	R = resistance;
	Ri = is_resistance;
	p0 = pressure_out;
	rho = density;

	// setting type_code based on type
	if(type == "junction")
	{
		type_code = 0;
	}
	else if(type == "periferia")
	{
		type_code = 1;
	}
	else if(type == "heart")
	{
		type_code = 2;
	}
}
   
//--------------------------------------------------------------
vector<double> moc_node::boundary_coefficients()
{
	double num   = -p0/ R * Ri;
	double denum = -1./ R * Ri;
	
	vector<double> out{num,denum};

	return out;
}

//--------------------------------------------------------------
void moc_node::boundary_variables(double p)
{
	double Q = (p-p0)/R * Ri;

	if(do_save_memory)
	{
		pressure.push_back(p);
		volume_flow_rate.push_back(Q);
	}
}
