#include "moc_edge.h"

using namespace std;

//--------------------------------------------------------------
moc_edge::moc_edge(string a_ID)
{
	ID = a_ID;
}

//--------------------------------------------------------------
moc_edge::~moc_edge(){}

//--------------------------------------------------------------
void moc_edge::print_input()
{
	printf("\n %8s, %8s, %12s, %8s, %8s, %6.3f, %6.3f, %6.3f, %6.3f, %6.3f, %3i, %6.3e, %6.3e, %6.3e, %6.3e, %6.3e", "vis", ID.c_str(), name.c_str(), node_name_start.c_str(), node_name_end.c_str(), dns, dne, sns, sne, l, nx, E1, E2, eta2, Rs, Re);
}

//--------------------------------------------------------------
void moc_edge::initialization(double pressure_initial)
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
   volume_flow_rate_start.clear();
   volume_flow_rate_end.clear();
   mass_flow_rate_start.clear();
   mass_flow_rate_end.clear();

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
	p.assign(nx,pressure_initial);
	v.assign(nx,1.);
	for(int i=0; i<nx; i++)
	{	
		// linear interpolation
		double nx_d = (double)nx;
		double i_d = (double)i;
		double dnp = dns*(nx_d-i_d-1.)/(nx_d-1.) + dne*i_d/(nx_d-1.);
		double snp = sns*(nx_d-i_d-1.)/(nx_d-1.) + sne*i_d/(nx_d-1.);

		a[i] = pow(E1*snp/(dnp*rho),0.5);
		d[i] = dnp;
		A[i] = dnp*dnp*pi*0.25;
	}
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
	if(do_save_memory)
	{
		save_field_variables();
	}
}

//--------------------------------------------------------------
void moc_edge::set_pressure_upstream(double p_in)
{
	p[0] = p_in;
}

//--------------------------------------------------------------
double moc_edge::new_timestep()
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
	for(unsigned int i=0; i<nx; i++) // TODO: i=1 tol eleg
	{
		if(dt[i] < dt_min)
		{
			dt_min = dt[i];
		}
		if(dt[i]<0)
		{
			cout << " Negative time step in edge: " << name << " , at i " << i << "-th inner node, dt: " << dt[i] << endl;
			cout <<  " velocity:      " << v[i] << endl;
			cout <<  " wave velocity: " << a[i] << endl;
			cout <<  " pressure:      " << p[i] << endl;
		}
	}

	return dt_min;
} 

//--------------------------------------------------------------
void moc_edge::solve()
{	
	for(unsigned int i=1; i<nx-1; i++)
	{
		xp[i] = x[i-1] + dt[i]*(v[i-1]+a[i-1]);

		// new field variables
		vp[i] = (p[i-1]-p[i+1] + rho*(a[i-1]*v[i-1] + a[i+1]*v[i+1]) - dt[i]*rho*(a[i-1]*JL(i) + a[i+1]*JR(i))) / (rho*(a[i-1]+a[i+1]));
		pp[i] = p[i+1] + rho*a[i+1]*(vp[i]-v[i+1]) + rho*a[i+1]*dt[i]*JR(i);

		if(vp[i]>.9*a[i])
		{
			vp[i] = .9*a[i];
		}
		else if(vp[i]<-.9*a[i])
		{
			vp[i] = -.9*a[i];
		}
	}
}

//--------------------------------------------------------------
void moc_edge::interpolate(double dt_real)
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
vector<vector<double> > moc_edge::backward_solver(vector<double> t_d, vector<double> p_d, vector<double> vfr_d)
{	
	// giving initial condition at the outlet in time
	initialization_back(t_d, p_d, vfr_d);

	while(x_back.back() < l)
	{	
		double dx_real = new_spacestep_back();
		solve_back();
		interpolate_back(dx_real);

		//update_variables(dt_back);
		double ex = exp(-E2/eta2*dt_back);
		for(unsigned int i=0; i<nt_back; i++)
		{
			if(i>0)
	      {
      		update_ith_variables(i,ex,p[i],epsz2[i-1],epsz[i-1]);
	      }
	      else
	      {
	      	update_ith_variables(i,ex,p[0],epsz2[0],epsz[0]);
	      }
		}
		reduce_field_vectors();
	}

	vector<vector<double> > out;
	out.push_back(t_back);
	out.push_back(p);

	return out;
}

//--------------------------------------------------------------
void moc_edge::initialization_back(vector<double> t_in, vector<double> p_in, vector<double> vfr_in)
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
   volume_flow_rate_start.clear();
   volume_flow_rate_end.clear();
   mass_flow_rate_start.clear();
   mass_flow_rate_end.clear();

   // matching the parameters for the short notations
	set_short_parameters();

	// backward calculation vars
	x_back.clear(); x_back.push_back(0.);
	dx_back.clear();

	// number of timesteps
	// nt_back = (int) ceil((t_in.back() - t_in[0])/dt_back_max);
	nt_back = ane*nx*(t_in.back() - t_in[0])/l;

	// time step
	dt_back = (t_in.back() - t_in[0]) / ((double)nt_back-1.);

	// interpolating to equidistant time coordinates
	t_back.clear(); t_back.resize(nt_back);
	for(unsigned int i=0; i<nt_back; i++)
	{
		t_back[i] = t_in[0] + dt_back * i;
	}

	// interpolating pressure and velocity
	p.clear();     p.resize(nt_back);
	v.clear();     v.resize(nt_back);
	epsz2.clear();	epsz2.resize(nt_back);
	epsz.clear();	epsz.resize(nt_back);
	a.clear();		a.resize(nt_back);
	d.clear();		d.resize(nt_back);
	A.clear();		A.resize(nt_back);

	// giving initial conditions
	a[0]		= pow(E1*sne/(dne*rho),0.5);
	d[0]		= dne;
	A[0]		= dne*dne*pi/4.;
	epsz[0]	= 0.;
	epsz2[0]	= 0.;

	int j=0;
	for(unsigned int i=0; i<nt_back; i++)
	{
		double ti= t_back[i];

      // finding the position for linear interpolation
      bool got_it = false;
      while(!got_it)
      {
         if(ti >= t_in[j] && ti <= t_in[j+1])
         {
            got_it=true;
         }
         else
         {
            j++;
         }
      }
      double p_h = p_in[j+1];
      double p_l = p_in[j];
      double vfr_h = vfr_in[j+1];
      double vfr_l = vfr_in[j];
      double t_h = t_in[j+1];
      double t_l = t_in[j];

      p[i] = (p_h-p_l)/(t_h-t_l) * (ti-t_l) + p_l;
      if(i>0)
      {
			double ex = exp(-E2/eta2*dt_back);
      	update_ith_variables(i,ex,p[i],epsz2[i-1],epsz[i-1]);
	      v[i] = ((vfr_h-vfr_l)/(t_h-t_l) * (ti-t_l) + vfr_l) / A[i];
      }
      else
      {
   	   v[i] = ((vfr_h-vfr_l)/(t_h-t_l) * (ti-t_l) + vfr_l) / Ane;
      }
	}

	// clearing everything and resizing
	dt.clear(); 	dt.resize(nt_back);
	pp.clear();    pp.resize(nt_back);
	pq.clear();    pq.resize(nt_back);
	vp.clear();    vp.resize(nt_back);
	vq.clear();    vq.resize(nt_back);
	x.clear();		
	xp.clear();		xp.resize(nt_back);
	// geodetic differences are cleard as x coordinates are not known
	h.clear();
}

//--------------------------------------------------------------
double moc_edge::new_spacestep_back()
{	
	// reducing the size of xp and x with 2
	nt_back -= 2;
	dx_back.clear(); dx_back.resize(nt_back);

	// inner points
	for(unsigned int i=1; i<nt_back+1; i++)
	{
		dx_back[i-1] = 0.5*dt_back * (-v[i-1]+a[i-1]+v[i+1]+a[i+1]);
	}

	// finding the minimum time step
	double dx_min = dx_back[0];
	int index=0;;
	for(unsigned int i=0; i<nt_back; i++)
	{
		if(dx_back[i] < dx_min)
		{
			dx_min = dx_back[i];
			index = i;
		}
	}

	if(dx_min<0.)
	{
      printf("\n !WARNING! space step is negative: %6.3e during BACKWARD calculation at moc_edge %-6s", dx_min, name.c_str());
      cout << endl;
	}

	// checking whether the solver overshoots in x
	double x_new = x_back.back() + dx_min;
	if(x_new >= l)
	{
		dx_min = l-x_back.back();
		x_new = l;
	}

	// saving x point
	x_back.push_back(x_new);

	return dx_min;
}

//--------------------------------------------------------------
void moc_edge::solve_back()
{	
	// resizing vectors
	pp.clear(); pp.resize(nt_back);
	vp.clear(); vp.resize(nt_back);
	tp_back.clear(); tp_back.resize(nt_back);

	double dt;
	for(unsigned int i=1; i<nt_back+1; i++)
	{
		// timstep
		dt = dx_back[i-1]/(-v[i-1]+a[i-1]);
		// new vars
		tp_back[i-1] = t_back[i-1] + dt;

		double x = x_back[x_back.size()-2];
		double h = hs + x/l*(he-hs);
		double xp = x + dx_back[i-1];
		double ja = JA(dt, p[i-1], v[i-1], a[i-1], epsz[i-1], epsz2[i-1], d[i-1], h, xp, x);
		double jb = JB(dt, p[i+1], v[i+1], a[i+1], epsz[i+1], epsz2[i+1], d[i+1], h, xp, x);

		pp[i-1] = rho*a[i-1]*a[i+1]/(a[i-1]+a[i+1]) * (dt*(ja+jb) + p[i-1]/(rho*a[i-1]) + p[i+1]/(rho*a[i+1]) - v[i-1] + v[i+1]);
		vp[i-1] = v[i-1] + (pp[i-1]-p[i-1])/(rho*a[i-1]) - dt*ja;
	}
}

//--------------------------------------------------------------
void moc_edge::interpolate_back(double dx_real)
{
	// resize vectors
	pq.clear(); pq.resize(nt_back);
	vq.clear(); vq.resize(nt_back);

	// interpolating in time
	for(unsigned int i=0; i<nt_back; i++)
	{
		if(tp_back[i]<t_back[i+1])
		{
			pq[i] = p[i] + (tp_back[i]-t_back[i])/(t_back[i+1]-t_back[i])*(p[i+1]-p[i]);
			vq[i] = v[i] + (tp_back[i]-t_back[i])/(t_back[i+1]-t_back[i])*(v[i+1]-v[i]);
		}
		else
		{
			pq[i] = p[i+1] + (tp_back[i]-t_back[i+1])/(t_back[i+2]-t_back[i+1])*(p[i+2]-p[i+1]);
			vq[i] = v[i+1] + (tp_back[i]-t_back[i+1])/(t_back[i+2]-t_back[i+1])*(v[i+2]-v[i+1]);
		}
	}	

	// interpolating in space
	for(unsigned int i=0; i<nt_back; i++)
	{
		pq[i] = pq[i] + dx_real/dx_back[i] * (pp[i]-pq[i]);
		vq[i] = vq[i] + dx_real/dx_back[i] * (vp[i]-vq[i]);
	}

	// resizing vectors
	v.clear(); v.resize(nt_back);
	p.clear(); p.resize(nt_back);

	// giving boundaries
	v[0] = vq[0]; v[nt_back-1] = vq[nt_back-1];
	p[0] = pq[0]; p[nt_back-1] = pq[nt_back-1];

	// new time step
	dt_back = (tp_back[nt_back-1]-tp_back[0]) / (nt_back-1);
	t_back.clear(); t_back.resize(nt_back);
	t_back[0] = tp_back[0];
	for(unsigned int i=1; i<nt_back; i++)
	{
		t_back[i] = t_back[i-1] + dt_back;
	}

	// interpolating back to equidistant time scale
	for(unsigned int i=1; i<nt_back-1; i++)
	{
		if(t_back[i]<tp_back[i])
		{
			v[i] = vq[i-1] + (t_back[i]-tp_back[i-1])/(tp_back[i]-tp_back[i-1])*(vq[i]-vq[i-1]);
			p[i] = pq[i-1] + (t_back[i]-tp_back[i-1])/(tp_back[i]-tp_back[i-1])*(pq[i]-pq[i-1]);
		}
		else
		{
			v[i] = vq[i] + (t_back[i]-tp_back[i])/(tp_back[i+1]-tp_back[i])*(vq[i+1]-vq[i]);
			p[i] = pq[i] + (t_back[i]-tp_back[i])/(tp_back[i+1]-tp_back[i])*(pq[i+1]-pq[i]);
		}
	}
}

//--------------------------------------------------------------
void moc_edge::reduce_field_vectors()
{
	epsz2.pop_back();
	epsz.pop_back();
	a.pop_back();
	d.pop_back();
	A.pop_back();
}

//--------------------------------------------------------------
double moc_edge::JL(int i)
{
	double out = JL(dt[i], p[i-1], v[i-1], a[i-1], epsz[i-1], epsz2[i-1], d[i-1], h[i-1], xp[i], x[i-1]);
	return out;
}

//--------------------------------------------------------------
double moc_edge::JL(double dt, double p, double v, double a, double epsz, double epsz2, double d, double h, double xp, double x)
{
	// temp variable
	double dnp = dns + xp/l*(dne-dns);
	double snp = sns + xp/l*(sne-sns);

	double AL = ( (p-p0)*dnp*(2.*epsz+1)/(eta2*2.*snp) - E2/eta2*epsz2 )*exp(-E2/eta2*dt);

	// geodetic height difference
	double hp = hs + xp/l*(he-hs);

	// output
	double JL = g*(hp-h)/(xp-x) + 32.*nu*v/(d*d) + 2*a/(2.*epsz+1.)*AL;

	return JL;
}

//--------------------------------------------------------------
double moc_edge::JR(int i)
{
	double out = JR(dt[i], p[i+1], v[i+1], a[i+1], epsz[i+1], epsz2[i+1], d[i+1], h[i+1], xp[i], x[i+1]);
	return out;
}

//--------------------------------------------------------------
double moc_edge::JR(double dt, double p, double v, double a, double epsz, double epsz2, double d, double h, double xp, double x)
{
	// temp variables
	double dnp = dns + xp/l*(dne-dns);
	double snp = sns + xp/l*(sne-sns);

	double AR = ( (p-p0)*dnp*(2.*epsz+1)/(eta2*2.*snp) - E2/eta2*epsz2 )*exp(-E2/eta2*dt);

	// geodetic height difference
	double hp = hs + xp/l*(he-hs);

	// output
	double JR = g*(hp-h)/(xp-x) + 32.*nu*v/(d*d) - 2.*a/(2.*epsz+1.)*AR;

	return JR;
}

//--------------------------------------------------------------
double moc_edge::JA(double dt, double p, double v, double a, double epsz, double epsz2, double d, double h, double xp, double x)
{
	return JR(dt, p, v, a, epsz, epsz2, d, h, xp, x);
}

//--------------------------------------------------------------
double moc_edge::JB(double dt, double p, double v, double a, double epsz, double epsz2, double d, double h, double xp, double x)
{
	return JL(dt, p, v, a, epsz, epsz2, d, h, xp, x);
}

//--------------------------------------------------------------
void moc_edge::update_variables(double dt)
{
	double ex = exp(-E2/eta2*dt);
	for(unsigned int i=0; i<epsz2.size(); i++)
	{	
		update_ith_variables(i,ex,p[i],epsz2[i],epsz[i]);
	}
}

//--------------------------------------------------------------
void moc_edge::update_ith_variables(int i, double ex, double p_new, double epsz2_old, double epsz_old)
{
	// linear interpolation of nominal diameter, thickness and sound velocity
	double dnp = dns + (dne-dns) * ((double)i/((double)nx-1.));
	double snp = sns + (sne-sns) * ((double)i/((double)nx-1.));
	double anp = sqrt(E1*snp / (rho*dnp));

	epsz2[i] = (p_new-p0)*(2.*epsz_old+1.)*(1.-ex) * dnp / (E2*2.*snp) + epsz2_old*ex;
	epsz[i]  = epsz2[i] + (p_new-p0)*(2.*epsz_old+1.) *dnp / (E1*2.*snp*pow(epsz_old+1.,beta));
	a[i]     = anp*pow(epsz[i]+1.,beta/2.);
	d[i]     = dnp*(epsz[i]+1.);
	A[i]     = d[i]*d[i]*pi/4.;
}

//--------------------------------------------------------------
void moc_edge::save_field_variables()
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

	double vf_s = v[0] * d[0]*d[0]*pi/4.;
	double vf_e = v[nx-1] * d[nx-1]*d[nx-1]*pi/4.;
	volume_flow_rate_start.push_back(vf_s);
	volume_flow_rate_end.push_back(vf_e);

	double mf_s = vf_s * rho;
	double mf_e = vf_e * rho;
	mass_flow_rate_start.push_back(mf_s);
	mass_flow_rate_end.push_back(mf_e);
}

//--------------------------------------------------------------
vector<double> moc_edge::boundary_master_start(double dt_master)
{	
	// location of R point
	double xR = boundary_start_position(dt_master);

	// interpolating fied variables to xR
	double aR = (a[1]-a[0]) / x[1] * xR + a[0];
	double vR = (v[1]-v[0]) / x[1] * xR + v[0];
	double pR = (p[1]-p[0]) / x[1] * xR + p[0];
	double dR = (d[1]-d[0]) / x[1] * xR + d[0];
	double AR = dR*dR*pi/4.;
	double hR = (h[1]-h[0]) / x[1] * xR + h[0];
	double epszR = (epsz[1]-epsz[0]) / x[1] * xR + epsz[0];
	double epsz2R = (epsz2[1]-epsz2[0]) / x[1] * xR + epsz2[0];

	// source term
	double J = JR(dt_master, pR, vR, aR, epszR, epsz2R, dR, hR, xR, x[0]);

	// coefficients aq*qp + ap*pp = b
	double aq = 1./A[0];
	double ap = -1./(rho*aR);
	double b  = vR + ap*pR - dt_master*J;

	vector<double> out{aq,ap,b};
	return out;
}

//--------------------------------------------------------------
vector<double> moc_edge::boundary_master_end(double dt_master)
{	
	// location of R point
	double xL = boundary_end_position(dt_master);

	// interpolating fied variables to xR
	double aL = (a[nx-2]-a[nx-1]) / (x[nx-2]-x[nx-1]) * (xL - x[nx-1]) + a[nx-1];
	double vL = (v[nx-2]-v[nx-1]) / (x[nx-2]-x[nx-1]) * (xL - x[nx-1]) + v[nx-1];
	double pL = (p[nx-2]-p[nx-1]) / (x[nx-2]-x[nx-1]) * (xL - x[nx-1]) + p[nx-1];
	double dL = (d[nx-2]-d[nx-1]) / (x[nx-2]-x[nx-1]) * (xL - x[nx-1]) + d[nx-1];
	double AL = dL*dL*pi/4.;
	double hL = (h[nx-2]-h[nx-1]) / (x[nx-2]-x[nx-1]) * (xL - x[nx-1]) + h[nx-1];
	double epszL = (epsz[nx-2]-epsz[nx-1]) / (x[nx-2]-x[nx-1]) * (xL - x[nx-1]) + epsz[nx-1];
	double epsz2L = (epsz2[nx-2]-epsz2[nx-1]) / (x[nx-2]-x[nx-1]) * (xL - x[nx-1]) + epsz2[nx-1];

	// source term
	double J = JL(dt_master, pL, vL, aL, epszL, epsz2L, dL, hL, xL, x[nx-1]);

	// coefficients aq*qp + ap*pp = b
	double aq = 1./A[nx-1];
	double ap = 1./(rho*aL);
	double b  = vL + ap*pL - dt_master*J;

	vector<double> out{aq,ap,b};
	return out;
}

//--------------------------------------------------------------
vector<double> moc_edge::boundary_start_coefficients(double dt)
{
	double xR = boundary_start_position(dt);

	// interpolating fied variables to xR
	double aR = (a[1]-a[0]) / x[1] * xR + a[0];
	double vR = (v[1]-v[0]) / x[1] * xR + v[0];
	double pR = (p[1]-p[0]) / x[1] * xR + p[0];
	double dR = (d[1]-d[0]) / x[1] * xR + d[0];
	double AR = dR*dR*pi/4.;
	double hR = (h[1]-h[0]) / x[1] * xR + h[0];
	double epszR = (epsz[1]-epsz[0]) / x[1] * xR + epsz[0];
	double epsz2R = (epsz2[1]-epsz2[0]) / x[1] * xR + epsz2[0];

	// source term
	double J = JR(dt, pR, vR, aR, epszR, epsz2R, dR, hR, xR, x[0]);

	// C- characteristics constants
	double C1m = 1. / (rho*aR + Rs*A[nx-1]) * A[0];
	double C2m = aR/(aR + Rs/rho*A[nx-1]) * (-dt*J + vR - pR/(rho*aR)) * A[0];

	vector<double> out{C1m,C2m};
	return out;
}

//--------------------------------------------------------------
vector<double> moc_edge::boundary_end_coefficients(double dt)
{
	double xL = boundary_end_position(dt);

	// interpolating fied variables to xR
	double aL = (a[nx-2]-a[nx-1]) / (x[nx-2]-x[nx-1]) * (xL - x[nx-1]) + a[nx-1];
	double vL = (v[nx-2]-v[nx-1]) / (x[nx-2]-x[nx-1]) * (xL - x[nx-1]) + v[nx-1];
	double pL = (p[nx-2]-p[nx-1]) / (x[nx-2]-x[nx-1]) * (xL - x[nx-1]) + p[nx-1];
	double dL = (d[nx-2]-d[nx-1]) / (x[nx-2]-x[nx-1]) * (xL - x[nx-1]) + d[nx-1];
	double AL = dL*dL*pi/4.;
	double hL = (h[nx-2]-h[nx-1]) / (x[nx-2]-x[nx-1]) * (xL - x[nx-1]) + h[nx-1];
	double epszL = (epsz[nx-2]-epsz[nx-1]) / (x[nx-2]-x[nx-1]) * (xL - x[nx-1]) + epsz[nx-1];
	double epsz2L = (epsz2[nx-2]-epsz2[nx-1]) / (x[nx-2]-x[nx-1]) * (xL - x[nx-1]) + epsz2[nx-1];

	// source term
	double J = JL(dt, pL, vL, aL, epszL, epsz2L, dL, hL, xL, x[nx-1]);

	// C+ characteristics
	double C1p = -1./(rho*aL + Re*A[nx-1]) * A[nx-1];
	double C2p = aL/(aL + Re/rho*A[nx-1])*(-dt*J + vL + pL/(rho*aL)) * A[nx-1];

	vector<double> out{C1p,C2p};
	return out;
}

//--------------------------------------------------------------
void moc_edge::boundary_start_variables(double dt, double p, double q)
{
	vp[0] = q/A[0];
	if(vp[0]>.9*a[0])
	{
		vp[0] = .9*a[0];
	}
	else if(vp[0]<-.9*a[0])
	{
		vp[0] = -.9*a[0];
	}
	pp[0] = p - Rs*q;
}

//--------------------------------------------------------------
void moc_edge::boundary_end_variables(double dt, double p, double q)
{
	vp[nx-1] = q/A[nx-1];
	if(vp[nx-1]>.9*a[nx-1])
	{
		vp[nx-1] = .9*a[nx-1];
	}
	else if(vp[nx-1]<-.9*a[nx-1])
	{
		vp[nx-1] = -.9*a[nx-1];
	}
	pp[nx-1] = p + Re*q;
}

//--------------------------------------------------------------
double moc_edge::boundary_periferia(double dt, double p_out)
{	
	double xL = boundary_end_position(dt);

	// interpolating fied variables to xR
	double aL = (a[nx-2]-a[nx-1]) / (x[nx-2]-x[nx-1]) * (xL - x[nx-1]) + a[nx-1];
	double vL = (v[nx-2]-v[nx-1]) / (x[nx-2]-x[nx-1]) * (xL - x[nx-1]) + v[nx-1];
	double pL = (p[nx-2]-p[nx-1]) / (x[nx-2]-x[nx-1]) * (xL - x[nx-1]) + p[nx-1];
	double dL = (d[nx-2]-d[nx-1]) / (x[nx-2]-x[nx-1]) * (xL - x[nx-1]) + d[nx-1];
	double AL = dL*dL*pi/4.;
	double hL = (h[nx-2]-h[nx-1]) / (x[nx-2]-x[nx-1]) * (xL - x[nx-1]) + h[nx-1];
	double epszL = (epsz[nx-2]-epsz[nx-1]) / (x[nx-2]-x[nx-1]) * (xL - x[nx-1]) + epsz[nx-1];
	double epsz2L = (epsz2[nx-2]-epsz2[nx-1]) / (x[nx-2]-x[nx-1]) * (xL - x[nx-1]) + epsz2[nx-1];

	// source term
	double J = JL(dt, pL, vL, aL, epszL, epsz2L, dL, hL, xL, x[nx-1]);

	double v_e;
	v_e = (vL - dt*J - (p_out-pL)/(rho*aL)) / (1.+Re/rho*A[nx-1]/aL);
	vp[nx-1] = v_e;

	double p_e;
	p_e = p_out + Re*A[nx-1]*v_e;
	pp[nx-1] = p_e;

	return v_e*A[nx-1];
}


//--------------------------------------------------------------
double moc_edge::boundary_start_position(double dt)
{
	double dv = (v[1]-v[0]) / x[1];
	double da = (a[1]-a[0]) / x[1];
	double xR = dt*(a[0]-v[0]) / (1+dt*dv-dt*da);
	return xR;
}

//--------------------------------------------------------------
double moc_edge::boundary_end_position(double dt)
{
	double dv = (v[nx-1]-v[nx-2]) / (x[nx-1]-x[nx-2]);
	double da = (a[nx-1]-a[nx-2]) / (x[nx-1]-x[nx-2]);
	double xL = (x[nx-1] + dt*dv*x[nx-2] - dt*v[nx-2] + dt*da*x[nx-2] - dt*a[nx-2]) / (1 + dt*dv + dt*da);
	return xL;
}

//--------------------------------------------------------------
double moc_edge::upstream_boundary_p(double dt, double p_in)
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
	double AR = dR*dR*pi/4.;

	// source term
	double J = JR(dt, pR, vR, aR, epszR, epsz2R, dR, hR, xR, x[0]);

	double v = (vR+1./(rho*aR)*(p_in-pR)-dt*J) / (1.+Rs/rho*AR/aR);
	vp[0] = v;

	double p = p_in - Rs*A[0]*v;
	pp[0] = p;

	return v*A[0];
}

//--------------------------------------------------------------
double moc_edge::upstream_boundary_q(double dt, double q_in)
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
	double AR = dR*dR*pi/4.;

	// source term
	double J = JR(dt, pR, vR, aR, epszR, epsz2R, dR, hR, xR, x[0]);

	// from BC
	double v = q_in/A[0];
	vp[0] = v;

	double p = pR + rho*aR*(dt*J + v - vR);
	pp[0] = p;

	return pp[0];
}

//--------------------------------------------------------------
void moc_edge::set_short_parameters()
{
   l    = length;
   dns   = nominal_diameter_start;
   dne   = nominal_diameter_end;
   sns   = nominal_thickness_start;
   sne   = nominal_thickness_end;
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

   Ans   = dns*dns*pi/4.;
   Ane   = dne*dne*pi/4.;

   ans   = sqrt(E1*sns / (rho*dns));
   ane   = sqrt(E1*sne / (rho*dne));
}

