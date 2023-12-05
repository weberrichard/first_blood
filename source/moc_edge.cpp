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
	printf("\n %8s, %8s, %12s, %8s, %8s, %6.3f, %6.3f, %6.3f, %6.3f, %6.3f, %3i, %6.3e, %6.3e, %6.3e", "vis", ID.c_str(), name.c_str(), node_name_start.c_str(), node_name_end.c_str(), Ans, Ane, sns, sne, l, nx, E, Rs, Re);
}

//--------------------------------------------------------------
void moc_edge::initialization(double pressure_initial, int mat_type)
{
	// clearing time variables
	pressure_start.clear();
   pressure_end.clear();
   velocity_start.clear();
   velocity_end.clear();
   wave_speed_start.clear();
   wave_speed_end.clear();
   area_start.clear();
   area_end.clear();
   volume_flow_rate_start.clear();
   volume_flow_rate_end.clear();
   mass_flow_rate_start.clear();
   mass_flow_rate_end.clear();
   time.clear();
   time.push_back(0.);

   // setting material properties
	material_type = mat_type;

   // matching the parameters for the short notations
	set_short_parameters();

	// clearing everything and resizing
	dt.clear();    dt.resize(nx);
	p.clear();     p.resize(nx);
	pnew.clear();  pnew.resize(nx);
	v.clear();     v.resize(nx);
	vnew.clear();  vnew.resize(nx);
	a.clear();     a.resize(nx);
	anew.clear();  anew.resize(nx);
	A.clear();     A.resize(nx);
	Anew.clear();  Anew.resize(nx);
	x.clear();     x.resize(nx);

	// calculating the space coordinates and dx
	dx = l/(nx-1);
	for(int i=0; i<nx; i++)
	{
		x[i] = i*dx;
	}

	// giving initial conditions
	p.assign(nx,pressure_initial);
	v.assign(nx,0.);
	for(int i=0; i<nx; i++)
	{	
		double d;
		A[i] = nominal_area(x[i],d);
		a[i] = wave_speed(x[i],A[i]);
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
bool moc_edge::new_timestep()
{
	// left boundary
	dt[0] = dx / (a[1]-v[1]);

	// inner points
	for(unsigned int i=1; i<nx-1; i++)
	{
		// dt[i] = 2.*dx / (v[i-1]+a[i-1]-v[i+1]+a[i+1]);
		dt[i] = min(dx/(a[i-1]+v[i-1]),dx/(a[i+1]-v[i+1]));
	}

	// right boundary
	dt[nx-1] = dx / (a[nx-2]+v[nx-2]);

	// finding the minimum time step
	double dt_min = dt[0];
	for(unsigned int i=0; i<nx; i++)
	{
		if(dt[i] < dt_min)
		{
			dt_min = dt[i];
		}
		if(dt[i]<0.)
		{
			cout << " Negative time step in edge: " << name << " , at i " << i << "-th node, dt: " << dt[i] << " time: " << time.back() << endl;
			cout <<  " velocity left:        " << v[i-1] << endl;
			cout <<  " velocity right:       " << v[i+1] << endl;
			cout <<  " wave velocity left:   " << a[i-1] << endl;
			cout <<  " wave velocity right:  " << a[i+1] << endl;
			cout <<  " pressure left:        " << p[i-1] << endl;
			cout <<  " pressure right:       " << p[i+1] << endl;
			return false;
		}
	}

	dt_act = cfl*dt_min;
	return true;
}

//--------------------------------------------------------------
void moc_edge::solve_maccormack()
{
	double dp_dA, dp_dx, dF_dx;
	double As = A[0] - dt_act/dx*(A[1]*v[1]-A[0]*v[0]);
	double vs = v[0] - dt_act/dx*(v[1]*v[1]*.5 + p[1]/rho - v[0]*v[0]*.5 - p[0]/rho) + dt_act*(-8.*pi*nu*nu_f/A[0]*v[0]); // TODO: add gravity
	double ps = pressure(x[0],As,dp_dA,dp_dx,dF_dx);

	double Fs_L_1 = As*vs;
	double Fs_L_2 = vs*vs*.5 + ps/rho;
	// McCormack scheme
	for(unsigned int i=1; i<nx-1; i++)
	{
		// predictor step Us = U_i - dt/dx*(F_i+1 - F_i) + dt*S_i
		As = A[i] - dt_act/dx*(A[i+1]*v[i+1]-A[i]*v[i]);
		vs = v[i] - dt_act/dx*(v[i+1]*v[i+1]*.5 + p[i+1]/rho - v[i]*v[i]*.5 - p[i]/rho) + dt_act*(-8.*pi*nu*nu_f/A[i]*v[i]); // TODO: add gravity to source
		ps = pressure(x[i],As,dp_dA,dp_dx,dF_dx);

		double Fs_i_1 = As*vs;
		double Fs_i_2 = vs*vs*.5 + ps/rho;
		// corrector step Un+1 = .5*(Us+U_n) + dt/(2dx)*(Fs_i - Fs_i-1) + dt/2*Ss_i
		Anew[i] = .5*(As+A[i]) - dt_act/(2.*dx)*(Fs_i_1-Fs_L_1);
		vnew[i] = .5*(vs+v[i]) - dt_act/(2.*dx)*(Fs_i_2-Fs_L_2) + dt_act*.5*(-8.*pi*nu*nu_f/As*vs);
		pnew[i] = pressure(x[i],Anew[i],dp_dA,dp_dx,dF_dx);
		anew[i] = wave_speed(x[i],Anew[i]);

		// refreshing left star variables
		Fs_L_1 = Fs_i_1;
		Fs_L_2 = Fs_i_2;
	}
}

//--------------------------------------------------------------
void moc_edge::solve_moc()
{
	for(int i=1; i<nx-1; i++)
	{
		double J1_L, J2_R;
		double W1_L = W1L(i, dt_act, J1_L);
		double W2_R = W2R(i, dt_act, J2_R);

		double sn_dx, An_dx, dp_dA, dp_dx, dF_dx;
		double sn = nominal_wall_thickness(x[i],sn_dx);
		double An = nominal_area(x[i],An_dx);

		vnew[i] = -J1_L*dt_act + W1_L - J2_R*dt_act + W2_R;
		if(material_type == 0)
		{
			Anew[i] = pow(.25*pow(2.*rho*An/(beta*sn),.5) * (-J1_L*dt_act + W1_L + J2_R*dt_act - W2_R) + pow(An,.25),4.);
		}
		else if(material_type == 1)
		{
			double dn = pow(An/pi,.5)*2.;
			double F = (k1*exp(k2*dn*.5)+k3)/(1.-nu_p*nu_p);
			Anew[i] = pow(.5*pow(rho/(2.*pow(An,.5)*F),.5) * (-W1_L + dt_act*J1_L + W2_R - dt_act*J2_R) + pow(An,-.25) ,-4);
		}
		pnew[i] = pressure(x[i],Anew[i],dp_dA,dp_dx,dF_dx);
		anew[i] = wave_speed(x[i],Anew[i]);
	}
}

//--------------------------------------------------------------
double moc_edge::JL(int i)
{
	double out = JL(dt[i], p[i-1], v[i-1], a[i-1], A[i-1], x[i-1]);
	return out;
}

//--------------------------------------------------------------
double moc_edge::JL(double dt, double p, double v, double a, double A, double xp)
{
	double dp_dA, dp_dx, dF_dx;
	pressure(xp,A,dp_dA,dp_dx,dF_dx);

	double sn_dx, An_dx;
	double sn = nominal_wall_thickness(xp,sn_dx);
	double An = nominal_area(xp,An_dx);

	double C1;
	if(material_type == 0)
	{
		C1 = -2.*pow(beta*.5/rho,.5)*((0.5*pow(1./(An*sn),.5)*sn_dx-0.5*pow(sn/pow(An,3.),.5)*An_dx)*(pow(A,0.25)-pow(An,0.25))-pow((sn/An),.5)*0.25*pow(An,-0.75)*An_dx);
	}
	else if(material_type == 1)
	{
		double dn = pow(An/pi,.5)*2.;
		double F  = (k1*exp(k2*dn*.5)+k3)/(1.-nu_p*nu_p);
		C1 = -pow(pow(An,.5)/(2.*rho*F),.5)*(pow(An,-.25)-pow(A,-.25))*dF_dx + pow(F/(8.*rho),.5)*pow(A,-.25)*pow(An,-.75)*An_dx;
	}
    
	double JL = nu_f*4.*pi*nu*v/A + .5*dp_dx/rho + (v+a)*C1; // TODO: add gravity
	// double JL = 4.*pi*nu*v/A + .5*dp_dx/rho; // TODO: add gravity

	return JL;
}

//--------------------------------------------------------------
double moc_edge::JR(int i)
{
	double out = JR(dt[i], p[i+1], v[i+1], a[i+1], A[i+1], x[i+1]);
	return out;
}

//--------------------------------------------------------------
double moc_edge::JR(double dt, double p, double v, double a, double A, double xp)
{
	double dp_dA, dp_dx, dF_dx;
	pressure(xp,A,dp_dA,dp_dx,dF_dx);

	double sn_dx, An_dx;
	double sn = nominal_wall_thickness(xp,sn_dx);
	double An = nominal_area(xp,An_dx);

	double C2;
	if(material_type == 0)
	{
		C2 = 2.*pow(beta*.5/rho,.5)*((0.5*pow(1./(An*sn),.5)*sn_dx-0.5*pow(sn/pow(An,3.),.5)*An_dx)*(pow(A,0.25)-pow(An,0.25)) - pow((sn/An),.5)*0.25*pow(An,-0.75)*An_dx);
	}
	else if(material_type == 1)
	{
		double dn = pow(An/pi,.5)*2.;
		double F  = (k1*exp(k2*dn*.5)+k3)/(1.-nu_p*nu_p);
		C2 = pow(pow(An,.5)/(2.*rho*F),.5)*(pow(An,-.25)-pow(A,-.25))*dF_dx - pow(F/(8.*rho),.5)*pow(A,-.25)*pow(An,-.75)*An_dx;
	}

	// double JR = 4.*pi*nu*v/A + .5*dp_dx/rho; // TODO: add gravity
	double JR = nu_f*4.*pi*nu*v/A + .5*dp_dx/rho + (v-a)*C2; // TODO: add gravity

	return JR;
}

//--------------------------------------------------------------
double moc_edge::nominal_area(double xp, double &An_dx)
{
	// linear
	An_dx = 1./l*(Ane-Ans); 
	double A = Ans + xp/l*(Ane-Ans);
	return A;
}

//--------------------------------------------------------------
double moc_edge::nominal_wall_thickness(double xp, double &sn_dx)
{
	// return pow(ss*ss + xp/l*(se*se-ss*ss),.5); // quadratic, i.e. linear in area
	
	// linear
	sn_dx = 1./l*(sne-sns);
	double s = sns + xp/l*(sne-sns);
	return s; 
}

//--------------------------------------------------------------
double moc_edge::wave_speed(double x, double A)
{
	double dp_dA, dp_dx, dF_dx;
	double p = pressure(x,A,dp_dA,dp_dx,dF_dx);
	return pow(A/rho*dp_dA,.5);
}

//--------------------------------------------------------------
double moc_edge::pressure(double x, double A, double &dp_dA, double &dp_dx, double &dF_dx)
{
	if(material_type == 0) // linear
	{
		double sn_dx, An_dx;
		double sn = nominal_wall_thickness(x,sn_dx);
		double An = nominal_area(x,An_dx);
		double sqrt_A = pow(A,.5);
		double sqrt_An = pow(An,.5);

		dp_dA = beta*sn/(An*2.*sqrt_A);
		dp_dx = beta/An*((sqrt_A-sqrt_An)*(sn_dx-sn/An*An_dx)-sn/(2.*sqrt_An)*An_dx);
		dF_dx = 0.;
		return beta*sn/An*(sqrt_A-sqrt_An) + p0;
	}
	else if(material_type == 1) // Olufsen
	{
		double sn_dx, An_dx;
		double sn = nominal_wall_thickness(x,sn_dx);
		double An = nominal_area(x,An_dx);
		double rn = pow(An/pi,.5);
		double rn_dx = 1./(2.*rn*pi)*An_dx;
		double sqrt_A = pow(A,.5);
		double sqrt_An = pow(An,.5);

		double F = (k1*exp(k2*rn)+k3)/(1.-nu_p*nu_p);
		dF_dx = 1./(1.-nu_p*nu_p)*k2*k1*exp(k2*rn)*rn_dx;
		dp_dA = F/2.*pow(A,-1.5)*sqrt_An;
		dp_dx = dF_dx*(1.-sqrt_An/sqrt_A) - F/(2.*sqrt_An*sqrt_A)*An_dx;
		return F*(1.-sqrt_An/sqrt_A) + p0;
	}
	else
	{
		cout << "Unknown material type = " << material_type << endl;
		return -1.;
	}
}

//--------------------------------------------------------------
double moc_edge::area(double x, double p)
{
	if(material_type == 0) // linear
	{
		double sn_dx, An_dx;
		double sn = nominal_wall_thickness(x,sn_dx);
		double An = nominal_area(x,An_dx);
		return pow((p-p0)*An/(sn*beta) + pow(An,.5),2.);
	}
	else if(material_type == 1) // olufsen
	{
		double sn_dx, An_dx;
		double sn = nominal_wall_thickness(x,sn_dx);
		double An = nominal_area(x,An_dx);
		double dn = pow(An/pi,.5)*2.;

		double F = (k1*exp(k2*dn*.5)+k3)/(1.-nu_p*nu_p);
		return An*pow(1.-(p-p0)/F,-2.);
	}
	else
	{
		cout << "Unknown material type = " << material_type << endl;
		return -1.;
	}
}

//--------------------------------------------------------------
double moc_edge::W1L(int j, double dt, double &J_L)
{
	// location of R point
	double xL = left_position(j,dt);

	double vL = (v[j-1]-v[j]) / (x[j-1]-x[j]) * (xL - x[j]) + v[j];
	double pL = (p[j-1]-p[j]) / (x[j-1]-x[j]) * (xL - x[j]) + p[j];
	double aL = (a[j-1]-a[j]) / (x[j-1]-x[j]) * (xL - x[j]) + a[j];
	double AL = (A[j-1]-A[j]) / (x[j-1]-x[j]) * (xL - x[j]) + A[j];

	J_L = JL(dt, pL, vL, aL, AL, xL);

	double d;
	double snL = nominal_wall_thickness(xL,d);
	double AnL = nominal_area(xL,d);

	double W1;
	if(material_type == 0)
	{
		W1 = .5*vL + 2.*pow(beta*snL/(2.*rho*AnL),.5)*(pow(AL,.25)-pow(AnL,.25));
	}
	else if(material_type == 1)
	{
		double rnL = pow(AnL/pi,.5);
		double F = (k1*exp(k2*rnL)+k3)/(1.-nu_p*nu_p);
		W1 = .5*vL + pow(2.*pow(AnL,.5)*F/rho,.5) * (pow(AnL,-.25)-pow(AL,-.25));
	}

	return W1;
}

//--------------------------------------------------------------
double moc_edge::W2R(int j, double dt, double &J_R)
{
	// location of R point
	double xR = right_position(j,dt);

	double vR = (v[j+1]-v[j]) / (x[j+1]-x[j]) * (xR - x[j]) + v[j];
	double pR = (p[j+1]-p[j]) / (x[j+1]-x[j]) * (xR - x[j]) + p[j];
	double aR = (a[j+1]-a[j]) / (x[j+1]-x[j]) * (xR - x[j]) + a[j];
	double AR = (A[j+1]-A[j]) / (x[j+1]-x[j]) * (xR - x[j]) + A[j];

	J_R = JR(dt, pR, vR, aR, AR, xR);

	double d;
	double snR = nominal_wall_thickness(xR,d);
	double AnR = nominal_area(xR,d);

	double W2;
	if(material_type == 0)
	{
		W2 = .5*vR - 2.*pow(beta*snR/(2.*rho*AnR),.5)*(pow(AR,.25)-pow(AnR,.25));
	}
	else if(material_type == 1)
	{
		double rnR = pow(AnR/pi,.5);
		double F = (k1*exp(k2*rnR)+k3)/(1.-nu_p*nu_p);
		W2 = .5*vR - pow(2.*pow(AnR,.5)*F/rho,.5) * (pow(AnR,-.25)-pow(AR,-.25));
	}
	return W2;
}

//--------------------------------------------------------------
double moc_edge::right_position(int j, double dt)
{
	double dv = (v[j+1]-v[j]) / (x[j+1]-x[j]);
	double da = (a[j+1]-a[j]) / (x[j+1]-x[j]);
	double xR = x[j] + dt*(a[j]-v[j]) / (1.+dt*dv-dt*da);
	return xR;
}

//--------------------------------------------------------------
double moc_edge::left_position(int j, double dt)
{
	double dv = (v[j]-v[j-1]) / (x[j]-x[j-1]);
	double da = (a[j]-a[j-1]) / (x[j]-x[j-1]);
	double xL = x[j] - dt*(v[j]+a[j]) / (1.+dt*da+dt*dv);
	return xL;
}

//--------------------------------------------------------------
void moc_edge::save_field_variables()
{
	velocity_start.push_back(v[0]);
	velocity_end.push_back(v[nx-1]);
	pressure_start.push_back(p[0]);
	pressure_end.push_back(p[nx-1]);
	wave_speed_start.push_back(a[0]);
	wave_speed_end.push_back(a[nx-1]);
	area_start.push_back(A[0]);
	area_end.push_back(A[nx-1]);

	double vf_s = v[0] * A[0];
	double vf_e = v[nx-1] * A[nx-1];
	volume_flow_rate_start.push_back(vf_s);
	volume_flow_rate_end.push_back(vf_e);

	// debug TODO: del
	/*if(time.back()>-1.)
	{
		for(int i=0; i<nx; i++)
			printf("%3i, vfr: %6.3e, a: %6.3f, d: %8.5f, p: %6.3f, v: %8.5f\n",i,v[i]*A[i],a[i],d[i],p[i],v[i]);
		cout << endl;
		cin.get();
	}*/

	double mf_s = vf_s * rho;
	double mf_e = vf_e * rho;
	mass_flow_rate_start.push_back(mf_s);
	mass_flow_rate_end.push_back(mf_e);
}

//--------------------------------------------------------------
void moc_edge::update()
{
	p = pnew;
	v = vnew;
	a = anew;
	A = Anew;
	time.push_back(time.back() + dt_act);

	/*cout << "dt: " << dt_act << endl;
	for(int i=0; i<nx; i++)
	{
		printf("i: %3i, A: %6.3e, v: %6.3f, p: %8.1f, a: %6.3f\n",i,Anew[i],vnew[i],pnew[i],anew[i]);
	}
	cout << endl;
	cin.get();*/
}

//--------------------------------------------------------------
void moc_edge::set_newton_size(int n1, int n2)
{
	// memory allocation
	Jac.push_back(MatrixXd::Zero(2+n1,2+n1));
	y.push_back(VectorXd::Zero(2+n1));
	f.push_back(VectorXd::Zero(2+2*n1));

	// memory allocation
	Jac.push_back(MatrixXd::Zero(2+n2,2+n2));
	y.push_back(VectorXd::Zero(2+n2));
	f.push_back(VectorXd::Zero(2+n2));
}

//--------------------------------------------------------------
double moc_edge::initialization_newton_start()
{
	return v[0]*A[0];
}

//--------------------------------------------------------------
double moc_edge::initialization_newton_end()
{
	return v[nx-1]*A[nx-1];
}

//--------------------------------------------------------------
void moc_edge::set_short_parameters()
{
   l      = length;
   dns    = nominal_diameter_start;
   dne    = nominal_diameter_end;
   Ans    = dns*dns*pi*.25;
   Ane    = dne*dne*pi*.25;
   sns    = nominal_thickness_start;
   sne    = nominal_thickness_end;
   E      = elasticity;
   Rs     = resistance_start;
   Re     = resistance_end;
   nx     = division_points;
   g      = gravity;
   rho    = density;
   nu     = kinematic_viscosity;
   nu_f   = kinematic_viscosity_factor;
   hs     = geodetic_height_start;
   he     = geodetic_height_end;
   p0     = atmospheric_pressure;
   nu_p   = poisson_coefficient;
   cfl    = courant_number;

   if(material_type == 1) // olufsen model
   {
   	k1 = material_const[0];
		k2 = material_const[1];
		k3 = material_const[2];
   }

   beta = pow(pi,.5)*E/(1.-nu_p*nu_p);
}

//--------------------------------------------------------------
void moc_edge::print_vars()
{
	printf(" time: %8.5f\n",time.back());
	printf(" pressure , velocity ,     a    \n");
	for(int i=0; i<nx; i++)
	{
		printf(" %8.3f, %8.3f , %8.3f \n",p[i],v[i],a[i]);
	}
	printf("\n");
}

//--------------------------------------------------------------
void moc_edge::save_initials(FILE* out_file)
{
   double ps = p[0];
   double pe = p[nx-1];
   double vs = v[0];
   double ve = v[nx-1];
   double As = A[0];
   double Ae = A[nx-1];
   double as = a[0];
   double ae = a[nx-1];

   fprintf(out_file, "%8s, %9.7e, %9.7e, %9.7e, %9.7e, %9.7e, %9.7e, %9.7e, %9.7e\n",ID.c_str(),ps,pe,vs,ve,As,Ae,as,ae);
}

//--------------------------------------------------------------
void moc_edge::set_initials(vector<double> ic)
{
	pressure_start.clear();
   pressure_end.clear();
   velocity_start.clear();
   velocity_end.clear();
   wave_speed_start.clear();
   wave_speed_end.clear();
   area_start.clear();
   area_end.clear();
   volume_flow_rate_start.clear();
   volume_flow_rate_end.clear();
   mass_flow_rate_start.clear();
   mass_flow_rate_end.clear();

	p.clear();      p.resize(nx);
	v.clear();      v.resize(nx);
	a.clear();      a.resize(nx);
	A.clear();      A.resize(nx);

	for(int i=0; i<nx; i++)
	{	
		double a0 = 1. - (double)i/((double)nx-1.);
		double a1 = (double)i / ((double)nx - 1.);

		p[i] = a0*ic[0] + a1*ic[1];
		v[i] = a0*ic[2] + a1*ic[3];
		A[i] = a0*ic[4] + a1*ic[5];
		a[i] = a0*ic[6] + a1*ic[7];
	}

	// saving initial conditions
	if(do_save_memory)
	{
		save_field_variables();
	}
}

//----------------------------------------------------------------\\
// *** BOUNDARIES *** BOUNDARIES *** BOUNDARIES *** BOUNDARIES*** \\
//----------------------------------------------------------------\\

//--------------------------------------------------------------
double moc_edge::boundary_start_position(double dt)
{
	double dv = (v[1]-v[0]) / x[1];
	double da = (a[1]-a[0]) / x[1];
	double xR = dt*(a[0]-v[0]) / (1.+dt*dv-dt*da);
	return xR;
}

//--------------------------------------------------------------
double moc_edge::boundary_end_position(double dt)
{
	double dv = (v[nx-1]-v[nx-2]) / (x[nx-1]-x[nx-2]);
	double da = (a[nx-1]-a[nx-2]) / (x[nx-1]-x[nx-2]);
	double xL = x[nx-1] - dt*(v[nx-1]+a[nx-1]) / (1.+dt*da+dt*dv);
	return xL;
}

//--------------------------------------------------------------
void moc_edge::boundary_substitute_start(double t_act, double p, double q)
{
	pnew[0] = p - Rs*q;
	Anew[0] = area(x[0],pnew[0]);
	anew[0] = wave_speed(x[0],Anew[0]);
	vnew[0] = q/Anew[0];
}

//--------------------------------------------------------------
void moc_edge::boundary_substitute_end(double t_act, double p, double q)
{
	pnew[nx-1] = p + Re*q;
	Anew[nx-1] = area(x[nx-1],pnew[nx-1]);
	anew[nx-1] = wave_speed(x[nx-1],Anew[nx-1]);
	vnew[nx-1] = q/Anew[nx-1];
}

//--------------------------------------------------------------
vector<double> moc_edge::boundary_junction_start(double pp, double t_act)
{
	// actual time step
	double dt = t_act - time.back();

	double J2_R;
	double W2_R = W2R(0, dt, J2_R);

	double Ap = area(x[0],pp);
	double dp = pp*0.001;
	double Ap2 = area(x[0],pp+dp);

	double q,q2;
	if(material_type == 0)
	{
		q   = Ap *(-2.*dt*J2_R + 2*W2_R + 4.*pow(beta*sns/(2.*rho*Ans),.5)*(pow(Ap ,.25)-pow(Ans,.25)));
		q2  = Ap2*(-2.*dt*J2_R + 2*W2_R + 4.*pow(beta*sns/(2.*rho*Ans),.5)*(pow(Ap2,.25)-pow(Ans,.25)));
	}
	else if(material_type == 1)
	{
		double F = (k1*exp(k2*dns*.5)+k3)/(1.-nu_p*nu_p);
		q  = Ap *(-2.*dt*J2_R + 2*W2_R + 2.*pow(2.*pow(Ans,.5)*F/rho,.5)*(pow(Ans,-.25)-pow(Ap ,-.25)));
		q2 = Ap2*(-2.*dt*J2_R + 2*W2_R + 2.*pow(2.*pow(Ans,.5)*F/rho,.5)*(pow(Ans,-.25)-pow(Ap2,-.25)));
	}

	double dq_dp = (q2-q)/dp;

	vector<double> out{q,dq_dp};
	return out;
}

//--------------------------------------------------------------
vector<double> moc_edge::boundary_junction_end(double pp, double t_act)
{
	// actual time step
	double dt = t_act - time.back();

	double J1_L;
	double W1_L = W1L(nx-1, dt, J1_L);

	double Ap = area(x[nx-1],pp);
	double dp = pp*0.001;
	double Ap2 = area(x[nx-1],pp+dp);

	double q,q2;
	if(material_type == 0)
	{
		q   = Ap *(-2.*dt*J1_L + 2*W1_L - 4.*pow(beta*sne/(2.*rho*Ane),.5)*(pow(Ap ,.25)-pow(Ane,.25)));
		q2  = Ap2*(-2.*dt*J1_L + 2*W1_L - 4.*pow(beta*sne/(2.*rho*Ane),.5)*(pow(Ap2,.25)-pow(Ane,.25)));
	}
	else if(material_type == 1)
	{
		double F = (k1*exp(k2*dne*.5)+k3)/(1.-nu_p*nu_p);
		q   = Ap *(-2.*dt*J1_L + 2*W1_L - 2.*pow(2.*pow(Ane,.5)*F/rho,.5)*(pow(Ane,-.25)-pow(Ap ,-.25)));
		q2  = Ap2*(-2.*dt*J1_L + 2*W1_L - 2.*pow(2.*pow(Ane,.5)*F/rho,.5)*(pow(Ane,-.25)-pow(Ap2,-.25)));
	}
	
	double dq_dp = (q2-q)/dp;

	vector<double> out{q,dq_dp};
	return out;
}

//--------------------------------------------------------------
vector<double> moc_edge::boundary_newton_start(double qp, double pp, double t_act)
{
	// actual time step
	double dt = t_act - time.back();

	double J2_R;
	double W2_R = W2R(0, dt, J2_R);

	double Ap = area(x[0],pp);
	double dp = pp*0.001;
	double Ap2 = area(x[0],pp+dp);

	// characteristic equation
	double f_char,f_char2;
	if(material_type == 0)
	{
		f_char  = .5*qp/Ap  - 2.*pow(beta*sns/(2.*rho*Ans),.5)*(pow(Ap ,.25)-pow(Ans,.25)) + dt*J2_R - W2_R;
		f_char2 = .5*qp/Ap2 - 2.*pow(beta*sns/(2.*rho*Ans),.5)*(pow(Ap2,.25)-pow(Ans,.25)) + dt*J2_R - W2_R;
	}
	else if(material_type == 1)
	{
		double F = (k1*exp(k2*dns*.5)+k3)/(1.-nu_p*nu_p);
		f_char  = .5*qp/Ap  - pow(2.*pow(Ans,.5)*F/rho,.5)*(pow(Ans,-.25)-pow(Ap ,-.25)) + dt*J2_R - W2_R;
		f_char2 = .5*qp/Ap2 - pow(2.*pow(Ans,.5)*F/rho,.5)*(pow(Ans,-.25)-pow(Ap2,-.25)) + dt*J2_R - W2_R;
	}

	// derivative to pp
	double dchar_dp = (f_char2-f_char)/dp;

	// derivative to qp
	double dchar_dq = .5/Ap;

	vector<double> out{f_char,dchar_dp,dchar_dq};
	return out;
}

//--------------------------------------------------------------
vector<double> moc_edge::boundary_newton_end(double qp, double pp, double t_act)
{
	// actual time step
	double dt = t_act - time.back();

	double J1_L;
	double W1_L = W1L(nx-1, dt, J1_L);

	double Ap = area(x[nx-1],pp);
	double dp = pp*0.001;
	double Ap2 = area(x[nx-1],pp+dp);

	// characteristic equation
	double f_char,f_char2;

	if(material_type == 0)
	{
		f_char  = .5*qp/Ap  + 2.*pow(beta*sne/(2.*rho*Ane),.5)*(pow(Ap ,.25)-pow(Ane,.25)) + dt*J1_L - W1_L;
		f_char2 = .5*qp/Ap2 + 2.*pow(beta*sne/(2.*rho*Ane),.5)*(pow(Ap2,.25)-pow(Ane,.25)) + dt*J1_L - W1_L;
	}
	else if(material_type == 1)
	{
		double F = (k1*exp(k2*dne*.5)+k3)/(1.-nu_p*nu_p);
		f_char  = .5*qp/Ap  + pow(2.*pow(Ane,.5)*F/rho,.5)*(pow(Ane,-.25)-pow(Ap ,-.25)) + dt*J1_L - W1_L;
		f_char2 = .5*qp/Ap2 + pow(2.*pow(Ane,.5)*F/rho,.5)*(pow(Ane,-.25)-pow(Ap2,-.25)) + dt*J1_L - W1_L;
	}

	// derivative to pp
	double dchar_dp = (f_char2-f_char)/dp;

	// derivative to qp
	double dchar_dq = .5/Ap;

	vector<double> out{f_char,dchar_dp,dchar_dq};
	return out;
}

//--------------------------------------------------------------
double moc_edge::boundary_pressure_start(double dt, double p_in)
{
	int k=0;
	double p_s = p[0], p_old = 0., v_s = 0.;
	double dp = p_s*0.00001;
	do
	{	
		p_old = p_s;
		double f1 = f_pressure_start(p_s   , p_in, dt, v_s);
		double f2 = f_pressure_start(p_s+dp, p_in, dt, v_s);
		double df_dp = (f2-f1)/dp;
		p_s = p_s - f1/df_dp;
		k++;
	}
	while(abs(p_s-p_old) > 1e-5 && k<100);

	if(k>=100)
	{
		cout << "Iteration did NOT converge at upstream_boundary_p, ID: " << ID << "\nExiting..." << endl;
		exit(-1);
	}

	// update v_s and d with final p_s
	pnew[0] = p_s;
	f_pressure_start(p_s, p_in, dt, v_s);
	vnew[0] = v_s;
	Anew[0] = area(x[0],pnew[0]);
	anew[0] = wave_speed(x[0],Anew[0]);

	double q_in = vnew[0]*Anew[0];

	return q_in;
}

//--------------------------------------------------------------
double moc_edge::f_pressure_start(double pp, double p_in, double dt, double &v_s)
{
	double J2_R;
	double W2_R = W2R(0, dt, J2_R);

	double Ap = area(x[0],pp);

	if(material_type == 0)
	{
		v_s = -2.*dt*J2_R + 2.*W2_R + 4.*pow(beta*sns/(2.*rho*Ans),.5)*(pow(Ap,.25)-pow(Ans,.25));
	}
	else if(material_type == 1)
	{
		double F = (k1*exp(k2*dns*.5)+k3)/(1.-nu_p*nu_p);
		v_s = -2.*dt*J2_R + 2.*W2_R + 2.*pow(2.*pow(Ans,.5)*F/rho,.5)*(pow(Ans,-.25)-pow(Ap,-.25));
	}

	double f = pp-p_in + Rs*Ap*v_s;

	return f;
}

//--------------------------------------------------------------
double moc_edge::boundary_pressure_end(double dt, double p_in)
{
	int k=0;
	double p_e = p[nx-1], p_old = 0., v_e = 0.;
	double dp = p_e*0.001;
	do
	{	
		p_old = p_e;
		double f1 = f_pressure_end(p_e   , p_in, dt, v_e);
		double f2 = f_pressure_end(p_e+dp, p_in, dt, v_e);
		double df_dp = (f2-f1)/dp;
		p_e = p_e - f1/df_dp;
		k++;
	}
	while((abs(p_e-p_old) > 1e-5) && k<100);

	if(k>=100)
	{
		cout << "Iteration did NOT converge at downstream_boundary_p, ID: " << ID << "\nExiting..." << endl;
		exit(-1);
	}

	// update v_s and d with final p_s
	pnew[nx-1] = p_e;
	f_pressure_end(p_e, p_in, dt, v_e);
	vnew[nx-1] = v_e;
	Anew[nx-1] = area(x[nx-1],pnew[nx-1]);
	anew[nx-1] = wave_speed(x[nx-1],Anew[nx-1]);

	double q_in = vnew[nx-1]*Anew[nx-1];

	return q_in;
}

//--------------------------------------------------------------
double moc_edge::f_pressure_end(double pp, double p_in, double dt, double &v_e)
{
	double J1_L;
	double W1_L = W1L(nx-1, dt, J1_L);

	double Ap = area(x[nx-1],pp);

	if(material_type == 0)
	{
		v_e = -2.*dt*J1_L + 2.*W1_L - 4.*pow(beta*sne/(2.*rho*Ane),.5)*(pow(Ap,.25)-pow(Ane,.25));
	}
	else if(material_type == 1)
	{
		double F = (k1*exp(k2*dne*.5)+k3)/(1.-nu_p*nu_p);
		v_e = -2.*dt*J1_L + 2.*W1_L - 2.*pow(2.*pow(Ane,.5)*F/rho,.5)*(pow(Ane,-.25)-pow(Ap,-.25));
	}

	double f = pp-p_in - Re*Ap*v_e;

	return f;
}

//--------------------------------------------------------------
double moc_edge::boundary_flowrate_start(double dt, double q_in)
{
	int k=0;
	double p_s = p[0], p_old = -10., v_s = 0.;
	double dp = p_s*0.001;
	do
	{	
		p_old = p_s;
		double f1 = f_flowrate_start(p_s   , q_in, dt, v_s);
		double f2 = f_flowrate_start(p_s+dp, q_in, dt, v_s);
		double df_dp = (f2-f1)/dp;
		p_s = p_s - f1/df_dp;
		k++;
	}
	while(abs(p_s-p_old) > 1e-5 && k<100);

	if(k>=100)
	{
		cout << "Iteration did NOT converge at upstream_boundary_q, ID: " << ID << "\nExiting..." << endl;
		exit(-1);
	}

	// update v_s and d with final p_s
	pnew[0] = p_s;
	f_flowrate_start(p_s, q_in, dt, v_s);
	vnew[0] = v_s;
	Anew[0] = area(x[0],pnew[0]);
	anew[0] = wave_speed(x[0],Anew[0]);

	return pnew[0];
}

//--------------------------------------------------------------
double moc_edge::f_flowrate_start(double pp, double q_in, double dt, double &v_s)
{
	double J2_R;
	double W2_R = W2R(0, dt, J2_R);

	double Ap = area(x[0],pp);

	v_s = q_in/Ap;

	double f;
	if(material_type == 0)
	{
		f = .5*v_s - 2.*pow(beta*sns/(2.*rho*Ans),.5)*(pow(Ap,.25)-pow(Ans,.25)) + dt*J2_R - W2_R;
	}
	else if(material_type == 1)
	{
		double F = (k1*exp(k2*dns*.5)+k3)/(1.-nu_p*nu_p);
		f = .5*v_s - pow(2.*pow(Ans,.5)*F/rho,.5)*(pow(Ans,-.25)-pow(Ap,-.25)) + dt*J2_R - W2_R;
	}

	return f;
}

//--------------------------------------------------------------
double moc_edge::boundary_flowrate_end(double dt, double q_in)
{
	int k=0;
	double p_e = p[nx-1], p_old = -10., v_e = 0.;
	double dp = p_e*0.001;
	do
	{	
		p_old = p_e;
		double f1 = f_flowrate_end(p_e   , q_in, dt, v_e);
		double f2 = f_flowrate_end(p_e+dp, q_in, dt, v_e);
		double df_dp = (f2-f1)/dp;
		p_e = p_e - f1/df_dp;
		k++;
	}
	while(abs(p_e-p_old) > 1e-5 && k<100);

	if(k>=100)
	{
		cout << "Iteration did NOT converge at downstream_boundary_p, ID: " << ID << "\nExiting..." << endl;
		exit(-1);
	}

	// update v_e with final p_e
	pnew[nx-1] = p_e;
	f_flowrate_end(p_e, q_in, dt, v_e);
	vnew[nx-1] = v_e;
	Anew[nx-1] = area(x[nx-1],pnew[nx-1]);
	anew[nx-1] = wave_speed(x[nx-1],Anew[nx-1]);

	return pnew[nx-1];
}

//--------------------------------------------------------------
double moc_edge::f_flowrate_end(double pp, double q_in, double dt, double &v_e)
{
	double J1_L;
	double W1_L = W1L(nx-1, dt, J1_L);

	double Ap = area(x[nx-1],pp);

	v_e = q_in/Ap;

	double f;
	if(material_type == 0)
	{
		f = .5*v_e + 2.*pow(beta*sne/(2.*rho*Ane),.5)*(pow(Ap,.25)-pow(Ane,.25)) + dt*J1_L - W1_L;
	}
	else if(material_type == 1)
	{
		double F = (k1*exp(k2*dne*.5)+k3)/(1.-nu_p*nu_p);
		f = .5*v_e + pow(2.*pow(Ane,.5)*F/rho,.5)*(pow(Ane,-.25)-pow(Ap,-.25)) + dt*J1_L - W1_L;
	}

	return f;
}

//--------------------------------------------------------------
double moc_edge::boundary_velocity_start(double dt, double v_in, double &q_in)
{
	double J2_R;
	double W2_R = W2R(0, dt, J2_R);

	vnew[0] = v_in;
	if(material_type == 0)
	{
		Anew[0] = pow(.5*pow(2.*rho*Ans/(beta*sns),.5)*(dt*J2_R - W2_R + .5*vnew[0]) + pow(Ans,.25),4.);
	}
	else if(material_type == 1)
	{
		double F = (k1*exp(k2*dns*.5)+k3)/(1.-nu_p*nu_p);
		Anew[0] = pow(pow(rho/(2.*pow(Ans,.5)*F),.5)*(-dt*J2_R + W2_R - .5*vnew[0]) + pow(Ans,-.25),-4.);
	}
	q_in = Anew[0]*vnew[0];
	double dp_dA, dp_dx, dF_dx;
	pnew[0] = pressure(x[0],Anew[0],dp_dA,dp_dx,dF_dx);
	anew[0] = wave_speed(x[0],Anew[0]);

	double p_in = Rs*q_in + pnew[0];

	return p_in;
}

//--------------------------------------------------------------
double moc_edge::boundary_velocity_end(double dt, double v_in, double &q_in)
{
	double J1_L;
	double W1_L = W1L(nx-1, dt, J1_L);

	vnew[nx-1] = v_in;
	if(material_type == 0)
	{
		Anew[nx-1] = pow(.5*pow(2.*rho*Ane/(beta*sne),.5)*(-dt*J1_L + W1_L - .5*vnew[nx-1]) + pow(Ane,.25),4.);    
	}
	else if(material_type == 1)
	{
		double F = (k1*exp(k2*dne*.5)+k3)/(1.-nu_p*nu_p);
		Anew[nx-1] = pow(pow(rho/(2.*pow(Ane,.5)*F),.5)*(dt*J1_L - W1_L + .5*vnew[nx-1]) + pow(Ane,-.25),-4.);
	}
	q_in = Anew[nx-1]*vnew[nx-1];
	double dp_dA, dp_dx, dF_dx;
	pnew[nx-1] = pressure(x[nx-1],Anew[nx-1],dp_dA,dp_dx,dF_dx);
	anew[nx-1] = wave_speed(x[nx-1],Anew[nx-1]);

	double p_in = -Re*q_in + pnew[nx-1];

	return p_in;
}
