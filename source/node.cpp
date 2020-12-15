#include "node.h"

using namespace std;

node::node(string a_name)
{
   name = a_name;
}

//--------------------------------------------------------------
node::~node(){}

//--------------------------------------------------------------
void node::initialization()
{
	// clearing time variables
   pressure.clear();
   volume_flow_rate.clear();

   // setting short versions of parameters
   set_short_parameters();

   // saving initial conditions
   pressure.push_back(pressure_out);
   volume_flow_rate.push_back(0.);
}

//--------------------------------------------------------------
void node::set_short_parameters()
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
vector<double> node::boundary_coefficients()
{
	double num   = -p0/ (R*rho) * Ri;
	double denum = -1./ (R*rho) * Ri;
	
	vector<double> out{num,denum};

	return out;
}

//--------------------------------------------------------------
void node::boundary_variables(double p)
{
	pressure.push_back(p);
	double Q = (p-p0)/(R*rho);
	volume_flow_rate.push_back(Q);
}
