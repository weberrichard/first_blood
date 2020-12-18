#include "edge.h"

using namespace std;

//--------------------------------------------------------------
edge::edge(string a_name)
{
	name = a_name;
}

//--------------------------------------------------------------
edge::~edge(){}

//--------------------------------------------------------------
void edge::print_input()
{
	printf("\n %8s, %8s, %8s, %8s, %6.3f, %6.3f, %6.3f, %3i, %6.3e, %6.3e, %6.3e, %6.3e, %6.3e", "vis", name.c_str(), node_name_start.c_str(), node_name_end.c_str(), dn, sn, l, nx, E1, E2, eta2, Rs, Re);
}


//--------------------------------------------------------------
void edge::initialization()
{
	// clearing time variables
	pressure_start.clear();
   pressure_end.clear();
   velocity_start.clear();
   velocity_end.clear();
   wave_velocity_start.clear();
   wave_velocity_end.clear();
   total_deformation_start.clear();
   total_deformation_end.clear();
   damper_deformation_start.clear();
   damper_deformation_end.clear();
   diameter_start.clear();
   diameter_end.clear();

   // matching the parameters for the short notations
	set_short_parameters();

	// clearing everything and resizing
	dt.clear();    dt.resize(nx);
	p.clear();     p.resize(nx);
	pp.clear();    pp.resize(nx);
	pq.clear();    pq.resize(nx);
	v.clear();     v.resize(nx);
	vp.clear();    vp.resize(nx);
	vq.clear();    vq.resize(nx);
	a.clear();     a.resize(nx);
	d.clear();     d.resize(nx);
	A.clear();     A.resize(nx);
	epsz.clear();  epsz.resize(nx);
	epsz2.clear(); epsz2.resize(nx);
	x.clear();     x.resize(nx);
	xp.clear();    xp.resize(nx);

	// giving initial conditions
	p.assign(nx,p0);
	v.assign(nx,0.);
	a.assign(nx,pow(E1*sn/(dn*rho),0.5));
	d.assign(nx,dn);
	A.assign(nx,dn*dn*M_PI/4.);
	epsz.assign(nx,0.);
	epsz2.assign(nx,0.);

	// calculating the space coordinates and dx
	dx = l/(nx-1);
	for(int i=0; i<nx; i++)
	{
		x[i] = i*dx;
	}

	// setting the geodetic height distribution
	h.resize(nx);
	for(unsigned int i=0; i<nx; i++)
	{
		h[i] = hs + x[i]/l * (he-hs);
	}

	// saving initial conditions
	save_field_variables();
}

//--------------------------------------------------------------
void edge::set_pressure_upstream(double p_in)
{
	p[0] = p_in;
}

//--------------------------------------------------------------
double edge::new_timestep()
{
	// inner points
	for(unsigned int i=1; i<nx-1; i++)
	{
		dt[i] = 2.*dx / (v[i-1]+a[i-1]-v[i+1]+a[i+1]);
	}

	// left boundary
	dt[0] = dx / (a[1]-v[1]);

	// right boundary
	dt[nx-1] = dx / (a[nx-1]+v[nx-1]);

	// finding the minimum time step
	double dt_min = dt[0];
	for(unsigned int i=1; i<nx-1; i++)
	{
		if(dt[i] < dt_min)
		{
			dt_min = dt[i];
		}
	}

	return dt_min;
} 

//--------------------------------------------------------------
void edge::solve()
{	
	for(unsigned int i=1; i<nx-1; i++)
	{
		xp[i] = x[i-1] + dt[i]*(v[i-1]+a[i-1]);

		// new field variables
		vp[i] = (p[i-1]-p[i+1] + rho*(a[i-1]*v[i-1] + a[i+1]*v[i+1]) - dt[i]*rho*(a[i-1]*JL(i) + a[i+1]*JR(i))) / (rho*(a[i-1]+a[i+1]));
		pp[i] = p[i+1] + rho*a[i+1]*(vp[i]-v[i+1]) + rho*a[i+1]*dt[i]*JR(i);
	}
}

//--------------------------------------------------------------
double edge::JL(int i)
{
	double out =	JL(dt[i], p[i-1], v[i-1], a[i-1], epsz[i-1], epsz2[i-1], d[i-1], h[i-1], xp[i], x[i-1]);
	return out;
}

//--------------------------------------------------------------
double edge::JL(double dt, double p, double v, double a, double epsz, double epsz2, double d, double h, double xp, double x)
{
	// temp variable
	double AL = ( (p-p0)*dn*(2.*epsz+1)/(eta2*2.*sn) - E2/eta2*epsz2 )*exp(-E2/eta2*dt);

	// geodetic height difference
	double hp = hs + xp/l*(he-hs);

	// output
	double JL = g*(hp-h)/(xp-x) + 32.*nu*v/(d*d) + 2*a/(2.*epsz+1.)*AL;

	return JL;
}

//--------------------------------------------------------------
double edge::JR(int i)
{
	double out = JR(dt[i], p[i+1], v[i+1], a[i+1], epsz[i+1], epsz2[i+1], d[i+1], h[i+1], xp[i], x[i+1]);
	return out;
}

//--------------------------------------------------------------
double edge::JR(double dt, double p, double v, double a, double epsz, double epsz2, double d, double h, double xp, double x)
{
	// temp variables
	double AR = ( (p-p0)*dn*(2.*epsz+1)/(eta2*2.*sn) - E2/eta2*epsz2 )*exp(-E2/eta2*dt);

	// geodetic height difference
	double hp = hs + xp/l*(he-hs);

	// output
	double JR = g*(hp-h)/(xp-x) + 32.*nu*v/(d*d) - 2*a/(2.*epsz+1.)*AR;

	return JR;
}

//--------------------------------------------------------------
void edge::interpolate(double dt_real)
{
	// setting back dt-s at boundaries
	dt[0]    = dt_real;
	dt[nx-1] = dt_real;

	// interpolation in space, position: xp, old time level
	vq[0]    = v[0];
	pq[0]    = p[0];
	vq[nx-1] = v[nx-1];
	pq[nx-1] = p[nx-1];
	for(unsigned int i=1; i<nx-1; i++)
	{
		if(xp[i]<x[i])
		{
			vq[i] = v[i-1] + (xp[i]-x[i-1])/dx*(v[i]-v[i-1]);
			pq[i] = p[i-1] + (xp[i]-x[i-1])/dx*(p[i]-p[i-1]);
		}
		else
		{
			vq[i] = v[i] + (xp[i]-x[i])/dx*(v[i+1]-v[i]);
			pq[i] = p[i] + (xp[i]-x[i])/dx*(p[i+1]-p[i]);
		}
	}

	// interpolation in time to new time level
	for(unsigned int i=0; i<nx; i++)
	{
		vq[i] = vq[i] + dt_real/dt[i] * (vp[i]-vq[i]);
		pq[i] = pq[i] + dt_real/dt[i] * (pp[i]-pq[i]);
	}

	// interpolation in space to equidistance grid
	v[0]    = vq[0];
	p[0]    = pq[0];
	v[nx-1] = vq[nx-1];
	p[nx-1] = pq[nx-1];
	for(unsigned int i=1; i<nx-1; i++)
	{
		if(x[i]<xp[i])
		{
			v[i] = vq[i-1] + (x[i]-xp[i-1])/(xp[i]-xp[i-1])*(vq[i]-vq[i-1]);
			p[i] = pq[i-1] + (x[i]-xp[i-1])/(xp[i]-xp[i-1])*(pq[i]-pq[i-1]);
		}
		else
		{
			v[i] = vq[i] + (x[i]-xp[i])/(xp[i+1]-xp[i])*(vq[i+1]-vq[i]);
			p[i] = pq[i] + (x[i]-xp[i])/(xp[i+1]-xp[i])*(pq[i+1]-pq[i]);
		}
	}
}

//--------------------------------------------------------------
void edge::update_field_variables(double dt)
{
	double ex = exp(-E2/eta2*dt);
	for(unsigned int i=0; i<nx; i++)
	{	
		epsz2[i] = (p[i]-p0)*(2.*epsz[i]+1.)*(1.-ex) * dn / (E2*2.*sn) + epsz2[i]*ex;
		epsz[i]  = epsz2[i] + (p[i]-p0)*(2.*epsz[i]+1.) *dn / (E1*2.*sn*pow(epsz[i]+1.,beta));
		a[i]     = an*pow(epsz[i]+1.,beta/2.);
		d[i]     = dn*(epsz[i]+1.);
		A[i]     = d[i]*d[i]*M_PI/4.;
	}
}

//--------------------------------------------------------------
void edge::save_field_variables()
{
	velocity_start.push_back(v[0]);
	velocity_end.push_back(v[nx-1]);
	pressure_start.push_back(p[0]);
	pressure_end.push_back(p[nx-1]);
	wave_velocity_start.push_back(a[0]);
	wave_velocity_end.push_back(a[nx-1]);
	total_deformation_start.push_back(epsz[0]);
	total_deformation_end.push_back(epsz[nx-1]);
	damper_deformation_start.push_back(epsz2[0]);
	damper_deformation_end.push_back(epsz2[nx-1]);
	diameter_start.push_back(d[0]);
	diameter_end.push_back(d[nx-1]);

	double vf_s = v[0] * d[0]*d[0]*M_PI/4.;
	double vf_e = v[nx-1] * d[nx-1]*d[nx-1]*M_PI/4.;
	volume_flow_rate_start.push_back(vf_s);
	volume_flow_rate_end.push_back(vf_e);

	double mf_s = vf_s * rho;
	double mf_e = vf_e * rho;
	mass_flow_rate_start.push_back(mf_s);
	mass_flow_rate_end.push_back(mf_e);
}

//--------------------------------------------------------------
vector<double> edge::boundary_start_coefficients(double dt)
{
	double xR = boundary_start_position(dt);

	// interpolating fied variables to xR
	double aR = (a[1]-a[0]) / x[1] * xR + a[0];
	double vR = (v[1]-v[0]) / x[1] * xR + v[0];
	double pR = (p[1]-p[0]) / x[1] * xR + p[0];
	double dR = (d[1]-d[0]) / x[1] * xR + d[0];
	double AR = dR*dR*M_PI/4.;
	double hR = (h[1]-h[0]) / x[1] * xR + h[0];
	double epszR = (epsz[1]-epsz[0]) / x[1] * xR + epsz[0];
	double epsz2R = (epsz2[1]-epsz2[0]) / x[1] * xR + epsz2[0];

	// source term
	double J = JR(dt, pR, vR, aR, epszR, epsz2R, dR, hR, xR, x[0]);

	// C- characteristics constants
	double C1m = 1. / (rho*aR + Rs*rho*A[nx-1]) * A[0];
	double C2m = aR/(aR + Rs*A[nx-1]) * (-dt*J + vR - pR/(rho*aR)) * A[0];

	vector<double> out{C1m,C2m};
	return out;
}

//--------------------------------------------------------------
vector<double> edge::boundary_end_coefficients(double dt)
{
	double xL = boundary_end_position(dt);

	// interpolating fied variables to xR
	double aL = (a[nx-2]-a[nx-1]) / (x[nx-2]-x[nx-1]) * (xL - x[nx-1]) + a[nx-1];
	double vL = (v[nx-2]-v[nx-1]) / (x[nx-2]-x[nx-1]) * (xL - x[nx-1]) + v[nx-1];
	double pL = (p[nx-2]-p[nx-1]) / (x[nx-2]-x[nx-1]) * (xL - x[nx-1]) + p[nx-1];
	double dL = (d[nx-2]-d[nx-1]) / (x[nx-2]-x[nx-1]) * (xL - x[nx-1]) + d[nx-1];
	double AL = dL*dL*M_PI/4.;
	double hL = (h[nx-2]-h[nx-1]) / (x[nx-2]-x[nx-1]) * (xL - x[nx-1]) + h[nx-1];
	double epszL = (epsz[nx-2]-epsz[nx-1]) / (x[nx-2]-x[nx-1]) * (xL - x[nx-1]) + epsz[nx-1];
	double epsz2L = (epsz2[nx-2]-epsz2[nx-1]) / (x[nx-2]-x[nx-1]) * (xL - x[nx-1]) + epsz2[nx-1];

	// source term
	double J = JL(dt, pL, vL, aL, epszL, epsz2L, dL, hL, xL, x[nx-1]);

	// C+ characteristics
	double C1p = -1./(rho*aL + Re*rho*A[nx-1]) * A[nx-1];
	double C2p = aL/(aL + Re*A[nx-1])*(-dt*J + vL + pL/(rho*aL)) * A[nx-1];

	vector<double> out{C1p,C2p};
	return out;
}

//--------------------------------------------------------------
void edge::boundary_start_variables(double dt, double p, double q)
{
	vp[0] = q/A[0];
	pp[0] = p - Rs*rho*q;
}

//--------------------------------------------------------------
void edge::boundary_end_variables(double dt, double p, double q)
{
	vp[nx-1] = q/A[nx-1];
	pp[nx-1] = p + Re*rho*q;
}

//--------------------------------------------------------------
double edge::boundary_periferia(double dt, double p_out)
{	
	double xL = boundary_end_position(dt);

	// interpolating fied variables to xR
	double aL = (a[nx-2]-a[nx-1]) / (x[nx-2]-x[nx-1]) * (xL - x[nx-1]) + a[nx-1];
	double vL = (v[nx-2]-v[nx-1]) / (x[nx-2]-x[nx-1]) * (xL - x[nx-1]) + v[nx-1];
	double pL = (p[nx-2]-p[nx-1]) / (x[nx-2]-x[nx-1]) * (xL - x[nx-1]) + p[nx-1];
	double dL = (d[nx-2]-d[nx-1]) / (x[nx-2]-x[nx-1]) * (xL - x[nx-1]) + d[nx-1];
	double AL = dL*dL*M_PI/4.;
	double hL = (h[nx-2]-h[nx-1]) / (x[nx-2]-x[nx-1]) * (xL - x[nx-1]) + h[nx-1];
	double epszL = (epsz[nx-2]-epsz[nx-1]) / (x[nx-2]-x[nx-1]) * (xL - x[nx-1]) + epsz[nx-1];
	double epsz2L = (epsz2[nx-2]-epsz2[nx-1]) / (x[nx-2]-x[nx-1]) * (xL - x[nx-1]) + epsz2[nx-1];

	// source term
	double J = JL(dt, pL, vL, aL, epszL, epsz2L, dL, hL, xL, x[nx-1]);

	double v_e;
	v_e = (vL - dt*J - (p_out-pL)/(rho*aL)) / (1.+Re*A[nx-1]/aL);
	vp[nx-1] = v_e;

	double p_e;
	p_e = p_out + Re*rho*A[nx-1]*v_e;
	pp[nx-1] = p_e;

	return v_e*A[nx-1];
}


//--------------------------------------------------------------
double edge::boundary_start_position(double dt)
{
	double dv = (v[1]-v[0]) / x[1];
	double da = (a[1]-a[0]) / x[1];
	double xR = dt*(a[0]-v[0]) / (1+dt*dv-dt*da);
	return xR;
}

//--------------------------------------------------------------
double edge::boundary_end_position(double dt)
{
	double dv = (v[nx-1]-v[nx-2]) / (x[nx-1]-x[nx-2]);
	double da = (a[nx-1]-a[nx-2]) / (x[nx-1]-x[nx-2]);
	double xL = (x[nx-1] + dt*dv*x[nx-2] - dt*v[nx-2] + dt*da*x[nx-2] - dt*a[nx-2]) / (1 + dt*dv + dt*da);
	return xL;
}

//--------------------------------------------------------------
double edge::upstream_boundary(double dt, double p_in)
{
	double xR = boundary_start_position(dt);

	// interpolating fied variables to xR
	double aR = (a[1]-a[0]) / x[1] * xR + a[0];
	double vR = (v[1]-v[0]) / x[1] * xR + v[0];
	double pR = (p[1]-p[0]) / x[1] * xR + p[0];
	double dR = (d[1]-d[0]) / x[1] * xR + d[0];
	double epszR = (epsz[1]-epsz[0]) / x[1] * xR + epsz[0];
	double epsz2R = (epsz2[1]-epsz2[0]) / x[1] * xR + epsz2[0];
	double hR = (h[1]-h[0]) / x[1] * xR + h[0];
	double AR = dR*dR*M_PI/4.;

	// source term
	double J = JR(dt, pR, vR, aR, epszR, epsz2R, dR, hR, xR, x[0]);

	double v = (vR+1./(rho*aR)*(p_in-pR)-dt*J) / (1.+Rs*AR/aR);
	vp[0] = v;

	double p = p_in - Rs*rho*A[0]*v;
	pp[0] = p;

	return v*A[0];
}

//--------------------------------------------------------------
void edge::set_short_parameters()
{
   l    = length;
   dn   = nominal_diameter;
   sn   = nominal_thickness;
   eta2 = viscosity;
   E1   = elasticity_spring;
   E2   = elasticity_voigt;
   Rs   = resistance_start;
   Re   = resistance_end;
   nx   = division_points;
   g    = gravity;
   rho  = density;
   nu   = kinematic_viscosity;
   hs   = geodetic_height_start;
   he   = geodetic_height_end;
   p0   = atmospheric_pressure;

   an   = sqrt(E1*sn / (rho*dn));
}

