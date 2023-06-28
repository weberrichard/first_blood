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
void moc_edge::initialization(double pressure_initial, int mat_type, vector<double> mat_const)
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

	// setting material properties
	material_type = mat_type;
	material_const = mat_const;
	beta = pow(pi,.5)*E/(1.-nu_p*nu_p);

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
			cout << " Negative time step in edge: " << name << " , at i " << i << "-th inner node, dt: " << dt[i] << endl;
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
	double dp_dA, dp_dx;
	double As = A[0] - dt_act/dx*(A[1]*v[1]-A[0]*v[0]);
	double vs = v[0] - dt_act/dx*(v[1]*v[1]*.5 + p[1]/rho - v[0]*v[0]*.5 - p[0]/rho) + dt_act*(-8.*pi*nu/A[0]*v[0]); // TODO: add gravity
	double ps = pressure(x[0],As,dp_dA,dp_dx);

	double Fs_L_1 = As*vs;
	double Fs_L_2 = vs*vs*.5 + ps/rho;
	// McCormack scheme
	for(unsigned int i=1; i<nx-1; i++)
	{
		// predictor step Us = U_i - dt/dx*(F_i+1 - F_i) + dt*S_i
		As = A[i] - dt_act/dx*(A[i+1]*v[i+1]-A[i]*v[i]);
		vs = v[i] - dt_act/dx*(v[i+1]*v[i+1]*.5 + p[i+1]/rho - v[i]*v[i]*.5 - p[i]/rho) + dt_act*(-8.*pi*nu/A[i]*v[i]); // TODO: add gravity to source
		ps = pressure(x[i],As,dp_dA,dp_dx);

		double Fs_i_1 = As*vs;
		double Fs_i_2 = vs*vs*.5 + ps/rho;
		// corrector step Un+1 = .5*(Us+U_n) + dt/(2dx)*(Fs_i - Fs_i-1) + dt/2*Ss_i
		Anew[i] = .5*(As+A[i]) - dt_act/(2.*dx)*(Fs_i_1-Fs_L_1);
		vnew[i] = .5*(vs+v[i]) - dt_act/(2.*dx)*(Fs_i_2-Fs_L_2) + dt_act*.5*(-8.*pi*nu/As*vs);
		pnew[i] = pressure(x[i],Anew[i],dp_dA,dp_dx);
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
		double W1_L = W1L(i, dt_act);
		double W2_R = W2R(i, dt_act);

		double J1_L = JL(i);
		double J2_R = JR(i);

		double sn_dx, An_dx, dp_dA, dp_dx;
		double sn = nominal_wall_thickness(x[i],sn_dx);
		double An = nominal_area(x[i],An_dx);

		vnew[i] = -J1_L*dt_act + W1_L - J2_R*dt_act + W2_R;
		Anew[i] = pow(.25*pow(2.*rho*An/(beta*sn),.5) * (-J1_L*dt_act + W1_L + J2_R*dt_act - W2_R) + pow(An,.25),4.);
		pnew[i] = pressure(x[i],Anew[i],dp_dA,dp_dx);
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
	double dp_dA, dp_dx;
	pressure(xp,A,dp_dA,dp_dx);

	double sn_dx, An_dx;
	double sn = nominal_wall_thickness(xp,sn_dx);
	double An = nominal_area(xp,An_dx);

	double C1 = -2.*pow(beta*.5/rho,.5)*((0.5*pow(1./(An*sn),.5)*sn_dx-0.5*pow(sn/pow(An,3.),.5)*An_dx)*(pow(A,0.25)-pow(An,0.25))-pow((sn/An),.5)*0.25*pow(An,-0.75)*An_dx);
    
	double JL = 4.*pi*nu*v/A + .5*dp_dx/rho + (v+a)*C1; // TODO: add gravity
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
	double dp_dA, dp_dx;
	pressure(xp,A,dp_dA,dp_dx);

	double sn_dx, An_dx;
	double sn = nominal_wall_thickness(xp,sn_dx);
	double An = nominal_area(xp,An_dx);

	double C2 = 2.*pow(beta*.5/rho,.5)*((0.5*pow(1./(An*sn),.5)*sn_dx-0.5*pow(sn/pow(An,3.),.5)*An_dx)*(pow(A,0.25)-pow(An,0.25))-pow((sn/An),.5)*0.25*pow(An,-0.75)*An_dx);

	// double JR = 4.*pi*nu*v/A + .5*dp_dx/rho; // TODO: add gravity
	double JR = 4.*pi*nu*v/A + .5*dp_dx/rho + (v-a)*C2; // TODO: add gravity

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
	double dp_dA, dp_dx;
	double p = pressure(x,A,dp_dA,dp_dx);
	return pow(A/rho*dp_dA,.5);
}

//--------------------------------------------------------------
double moc_edge::pressure(double x, double A, double &dp_dA, double &dp_dx)
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
		return beta*sn/An*(sqrt_A-sqrt_An) + p0;
	}
	else if(material_type == 1) // Olufsen
	{
		double sn_dx, An_dx;
		double sn = nominal_wall_thickness(x,sn_dx);
		double An = nominal_area(x,An_dx);
		double rn = pow(An/pi,.5);
		double rn_dx = 1./(2.*pow(pi*An,.5))*An_dx;
		double sqrt_A = pow(A,.5);
		double sqrt_An = pow(An,.5);

		double k1 = material_const[0];
		double k2 = material_const[1];
		double k3 = material_const[2];

		double F = 1./(1.-nu_p*nu_p)*(k1*exp(k2*rn)+k3);
		double dF_dx = 1./(1.-nu_p*nu_p)*k2*rn*k1*exp(k2*rn)*rn_dx;
		dp_dA = F/2.*pow(A,-1.5);
		dp_dx = dF_dx*(1.-sqrt_An/sqrt_A) - F/(2.*sqrt_An*sqrt_A) + F/2.*sqrt_An*pow(A,-1.5);
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
		double rn = pow(An/pi,.5);

		double k1 = material_const[0];
		double k2 = material_const[1];
		double k3 = material_const[2];

		return An*pow(1.-(p-p0)*(1.-nu_p*nu_p)/(k1*exp(k2*rn)+k3),-2.);
	}
	else
	{
		cout << "Unknown material type = " << material_type << endl;
		return -1.;
	}
}

//--------------------------------------------------------------
double moc_edge::W1L(int j, double dt)
{
	// location of R point
	double xL = left_position(j,dt);

	double vL = (v[j-1]-v[j]) / (x[j-1]-x[j]) * (xL - x[j]) + v[j];
	double pL = (p[j-1]-p[j]) / (x[j-1]-x[j]) * (xL - x[j]) + p[j];
	double aL = (a[j-1]-a[j]) / (x[j-1]-x[j]) * (xL - x[j]) + a[j];
	double AL = (A[j-1]-A[j]) / (x[j-1]-x[j]) * (xL - x[j]) + A[j];

	double J = JL(dt, pL, vL, aL, AL, xL);

	double d;
	double snL = nominal_wall_thickness(xL,d);
	double AnL = nominal_area(xL,d);

	double W1 = .5*vL + 2.*pow(beta*snL/(2.*rho*AnL),.5)*(pow(AL,.25)-pow(AnL,.25));

	return W1;
}

//--------------------------------------------------------------
double moc_edge::W2R(int j, double dt)
{
	// location of R point
	double xR = right_position(j,dt);

	double vR = (v[j+1]-v[j]) / (x[j+1]-x[j]) * (xR - x[j]) + v[j];
	double pR = (p[j+1]-p[j]) / (x[j+1]-x[j]) * (xR - x[j]) + p[j];
	double aR = (a[j+1]-a[j]) / (x[j+1]-x[j]) * (xR - x[j]) + a[j];
	double AR = (A[j+1]-A[j]) / (x[j+1]-x[j]) * (xR - x[j]) + A[j];

	double J = JR(dt, pR, vR, aR, AR, xR);

	double d;
	double snR = nominal_wall_thickness(xR,d);
	double AnR = nominal_area(xR,d);

	double W2 = .5*vR - 2.*pow(beta*snR/(2.*rho*AnR),.5)*(pow(AR,.25)-pow(AnR,.25));

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
			printf("%3i, vfr: %6.3e, a: %6.3f, d: %8.5f, p: %6.3f, v: %8.5f\n",i,v[i]*d[i]*d[i]*pi*.25,a[i],d[i],p[i],v[i]);
		cout << endl;
		//cin.get();
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

	//cout << "dt: " << dt_act << endl;
	//for(int i=0; i<nx; i++)
	//{
	//	printf("i: %3i, A: %6.3e, v: %6.3f, p: %8.1f, a: %6.3f\n",i,Anew[i],vnew[i],pnew[i],anew[i]);
	//}
	//cout << endl;
	//cin.get();
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

//---------------------------------------------------------------\\
// *** BOUNDARIES *** BOUNDARIES *** BOUNDARIES *** BOUNDARIES***\\
//---------------------------------------------------------------\\

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
	double dt = t_act - time.back();

	pnew[0] = p - Rs*q;
	Anew[0] = area(x[0],pnew[0]);
	anew[0] = wave_speed(x[0],Anew[0]);
	vnew[0] = q/Anew[0];
}

//--------------------------------------------------------------
void moc_edge::boundary_substitute_end(double t_act, double p, double q)
{
	double dt = t_act - time.back();

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

	// location of R point
	double xR = boundary_start_position(dt);

	double vR = (v[1]-v[0]) / x[1] * xR + v[0];
	double pR = (p[1]-p[0]) / x[1] * xR + p[0];
	double aR = (a[1]-a[0]) / x[1] * xR + a[0];
	double AR = (A[1]-A[0]) / x[1] * xR + A[0];

	double J = JR(dt, pR, vR, aR, AR, xR);

	double d;
	double snR = nominal_wall_thickness(xR,d);
	double AnR = nominal_area(xR,d);

	double W2R = .5*vR - 2.*pow(beta*snR/(2.*rho*AnR),.5)*(pow(AR,.25)-pow(AnR,.25));
	double Ap = area(x[0],pp);

	double q  = Ap*(-2.*dt*J + 2*W2R + 4.*pow(beta*sns/(2.*rho*Ans),.5)*(pow(Ap,.25)-pow(Ans,.25)));

	double dp = pp*0.001;
	double Ap2 = area(x[0],pp+dp);
	double q2  = Ap2*(-2.*dt*J + 2*W2R + 4.*pow(beta*sns/(2.*rho*Ans),.5)*(pow(Ap2,.25)-pow(Ans,.25)));
	double dq_dp = (q2-q)/dp;

	vector<double> out{q,dq_dp};
	return out;
}

//--------------------------------------------------------------
vector<double> moc_edge::boundary_junction_end(double pp, double t_act)
{
	// actual time step
	double dt = t_act - time.back();

	// location of R point
	double xL = boundary_end_position(dt);

	double vL = (v[nx-2]-v[nx-1]) / (x[nx-2]-x[nx-1]) * (xL - x[nx-1]) + v[nx-1];
	double pL = (p[nx-2]-p[nx-1]) / (x[nx-2]-x[nx-1]) * (xL - x[nx-1]) + p[nx-1];
	double aL = (a[nx-2]-a[nx-1]) / (x[nx-2]-x[nx-1]) * (xL - x[nx-1]) + a[nx-1];
	double AL = (A[nx-2]-A[nx-1]) / (x[nx-2]-x[nx-1]) * (xL - x[nx-1]) + A[nx-1];

	double J = JL(dt, pL, vL, aL, AL, xL);

	double d;
	double snL = nominal_wall_thickness(xL,d);
	double AnL = nominal_area(xL,d);

	double W1L = .5*vL + 2.*pow(beta*snL/(2.*rho*AnL),.5)*(pow(AL,.25)-pow(AnL,.25));
	double Ap = area(x[nx-1],pp);

	double q  = Ap*(-2.*dt*J + 2*W1L - 4.*pow(beta*sne/(2.*rho*Ane),.5)*(pow(Ap,.25)-pow(Ane,.25)));
	
	double dp = pp*0.001;
	double Ap2 = area(x[nx-1],pp+dp);
	double q2  = Ap2*(-2.*dt*J + 2*W1L - 4.*pow(beta*sne/(2.*rho*Ane),.5)*(pow(Ap2,.25)-pow(Ane,.25)));
	double dq_dp = (q2-q)/dp;

	vector<double> out{q,dq_dp};
	return out;
}

//--------------------------------------------------------------
vector<double> moc_edge::boundary_newton_start(double qp, double pp, double t_act)
{
	// actual time step
	double dt = t_act - time.back();

	// location of R point
	double xR = boundary_start_position(dt);

	double vR = (v[1]-v[0]) / x[1] * xR + v[0];
	double pR = (p[1]-p[0]) / x[1] * xR + p[0];
	double aR = (a[1]-a[0]) / x[1] * xR + a[0];
	double AR = (A[1]-A[0]) / x[1] * xR + A[0];

	double J = JR(dt, pR, vR, aR, AR, xR);

	double d;
	double snR = nominal_wall_thickness(xR,d);
	double AnR = nominal_area(xR,d);

	double W2R = .5*vR - 2.*pow(beta*snR/(2.*rho*AnR),.5)*(pow(AR,.25)-pow(AnR,.25));
	double Ap = area(x[0],pp);

	// characteristic equation
	double f_char = .5*qp/Ap - 2.*pow(beta*sns/(2.*rho*Ans),.5)*(pow(Ap,.25)-pow(Ans,.25)) + dt*J - W2R;

	// derivative to pp
	double dp = pp*0.001;
	double Ap2 = area(x[0],pp+dp);
	double f_char2 = .5*qp/Ap2 - 2.*pow(beta*sns/(2.*rho*Ans),.5)*(pow(Ap2,.25)-pow(Ans,.25)) + dt*J - W2R;
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

	// location of R point
	double xL = boundary_end_position(dt);

	double vL = (v[nx-2]-v[nx-1]) / (x[nx-2]-x[nx-1]) * (xL - x[nx-1]) + v[nx-1];
	double pL = (p[nx-2]-p[nx-1]) / (x[nx-2]-x[nx-1]) * (xL - x[nx-1]) + p[nx-1];
	double aL = (a[nx-2]-a[nx-1]) / (x[nx-2]-x[nx-1]) * (xL - x[nx-1]) + a[nx-1];
	double AL = (A[nx-2]-A[nx-1]) / (x[nx-2]-x[nx-1]) * (xL - x[nx-1]) + A[nx-1];

	double J = JL(dt, pL, vL, aL, AL, xL);

	double d;
	double snL = nominal_wall_thickness(xL,d);
	double AnL = nominal_area(xL,d);

	double W1L = .5*vL + 2.*pow(beta*snL/(2.*rho*AnL),.5)*(pow(AL,.25)-pow(AnL,.25));
	double Ap = area(x[nx-1],pp);

	// characteristic equation
	double f_char = .5*qp/Ap + 2.*pow(beta*sne/(2.*rho*Ane),.5)*(pow(Ap,.25)-pow(Ane,.25)) + dt*J - W1L;

	// derivative to pp
	double dp = pp*0.001;
	double Ap2 = area(x[nx-1],pp+dp);
	double f_char2 = .5*qp/Ap2 + 2.*pow(beta*sne/(2.*rho*Ane),.5)*(pow(Ap2,.25)-pow(Ane,.25)) + dt*J - W1L;
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
	double dp = p_s*0.001;
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
	double xR = boundary_start_position(dt);

	double vR = (v[1]-v[0]) / x[1] * xR + v[0];
	double pR = (p[1]-p[0]) / x[1] * xR + p[0];
	double aR = (a[1]-a[0]) / x[1] * xR + a[0];
	double AR = (A[1]-A[0]) / x[1] * xR + A[0];

	double J = JR(dt, pR, vR, aR, AR, xR);

	double d;
	double snR = nominal_wall_thickness(xR,d);
	double AnR = nominal_area(xR,d);

	double W2R = .5*vR - 2.*pow(beta*snR/(2.*rho*AnR),.5)*(pow(AR,.25)-pow(AnR,.25));
	double Ap = area(x[0],pp);

	v_s = -2.*dt_act*J + 2.*W2R + 4.*pow(beta*sns/(2.*rho*Ans),.5)*(pow(Ap,.25)-pow(Ans,.25));

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
	double xL = boundary_end_position(dt);

	double vL = (v[nx-2]-v[nx-1]) / (x[nx-2]-x[nx-1]) * (xL - x[nx-1]) + v[nx-1];
	double pL = (p[nx-2]-p[nx-1]) / (x[nx-2]-x[nx-1]) * (xL - x[nx-1]) + p[nx-1];
	double aL = (a[nx-2]-a[nx-1]) / (x[nx-2]-x[nx-1]) * (xL - x[nx-1]) + a[nx-1];
	double AL = (A[nx-2]-A[nx-1]) / (x[nx-2]-x[nx-1]) * (xL - x[nx-1]) + A[nx-1];

	double J = JL(dt, pL, vL, aL, AL, xL);

	double d;
	double snL = nominal_wall_thickness(xL,d);
	double AnL = nominal_area(xL,d);

	double W1L = .5*vL + 2.*pow(beta*snL/(2.*rho*AnL),.5)*(pow(AL,.25)-pow(AnL,.25));
	double Ap = area(x[nx-1],pp);

	v_e = -2.*dt_act*J + 2.*W1L - 4.*pow(beta*sns/(2.*rho*Ane),.5)*(pow(Ap,.25)-pow(Ane,.25));

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
		double xR = boundary_start_position(dt);

	double vR = (v[1]-v[0]) / x[1] * xR + v[0];
	double pR = (p[1]-p[0]) / x[1] * xR + p[0];
	double aR = (a[1]-a[0]) / x[1] * xR + a[0];
	double AR = (A[1]-A[0]) / x[1] * xR + A[0];

	double J = JR(dt, pR, vR, aR, AR, xR);

	double d;
	double snR = nominal_wall_thickness(xR,d);
	double AnR = nominal_area(xR,d);

	double W2R = .5*vR - 2.*pow(beta*snR/(2.*rho*AnR),.5)*(pow(AR,.25)-pow(AnR,.25));
	double Ap = area(x[0],pp);

	v_s = q_in/Ap;

	double f = .5*v_s - 2.*pow(beta*sns/(2.*rho*Ans),.5)*(pow(Ap,.25)-pow(Ans,.25)) + dt_act*J - W2R;

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
		double xL = boundary_end_position(dt);

	double vL = (v[nx-2]-v[nx-1]) / (x[nx-2]-x[nx-1]) * (xL - x[nx-1]) + v[nx-1];
	double pL = (p[nx-2]-p[nx-1]) / (x[nx-2]-x[nx-1]) * (xL - x[nx-1]) + p[nx-1];
	double aL = (a[nx-2]-a[nx-1]) / (x[nx-2]-x[nx-1]) * (xL - x[nx-1]) + a[nx-1];
	double AL = (A[nx-2]-A[nx-1]) / (x[nx-2]-x[nx-1]) * (xL - x[nx-1]) + A[nx-1];

	double J = JL(dt, pL, vL, aL, AL, xL);

	double d;
	double snL = nominal_wall_thickness(xL,d);
	double AnL = nominal_area(xL,d);

	double W1L = .5*vL + 2.*pow(beta*snL/(2.*rho*AnL),.5)*(pow(AL,.25)-pow(AnL,.25));
	double Ap = area(x[nx-1],pp);

	v_e = q_in/Ap;

	double f = .5*v_e + 2.*pow(beta*sne/(2.*rho*Ane),.5)*(pow(Ap,.25)-pow(Ane,.25)) + dt_act*J - W1L;

	return f;
}

//--------------------------------------------------------------
double moc_edge::boundary_velocity_start(double dt, double v_in, double &q_in)
{
	double xR = boundary_start_position(dt);

	double vR = (v[1]-v[0]) / x[1] * xR + v[0];
	double pR = (p[1]-p[0]) / x[1] * xR + p[0];
	double aR = (a[1]-a[0]) / x[1] * xR + a[0];
	double AR = (A[1]-A[0]) / x[1] * xR + A[0];

	double J = JR(dt, pR, vR, aR, AR, xR);

	double d;
	double snR = nominal_wall_thickness(xR,d);
	double AnR = nominal_area(xR,d);

	double W2R = .5*vR - 2.*pow(beta*snR/(2.*rho*AnR),.5)*(pow(AR,.25)-pow(AnR,.25));

	vnew[0] = v_in;
	Anew[0] = pow(.5*pow(2.*rho*Ans/(beta*sns),.5)*(dt_act*J - W2R + .5*vnew[0]) + pow(Ans,.25),4.);

	q_in = Anew[0]*vnew[0];
	double dp_dA, dp_dx;
	pnew[0] = pressure(x[0],Anew[0],dp_dA,dp_dx);
	anew[0] = wave_speed(x[0],Anew[0]);

	double p_in = Rs*q_in + pnew[0];

	return p_in;
}

//--------------------------------------------------------------
double moc_edge::boundary_velocity_end(double dt, double v_in, double &q_in)
{
	double xL = boundary_end_position(dt);

	double vL = (v[nx-2]-v[nx-1]) / (x[nx-2]-x[nx-1]) * (xL - x[nx-1]) + v[nx-1];
	double pL = (p[nx-2]-p[nx-1]) / (x[nx-2]-x[nx-1]) * (xL - x[nx-1]) + p[nx-1];
	double aL = (a[nx-2]-a[nx-1]) / (x[nx-2]-x[nx-1]) * (xL - x[nx-1]) + a[nx-1];
	double AL = (A[nx-2]-A[nx-1]) / (x[nx-2]-x[nx-1]) * (xL - x[nx-1]) + A[nx-1];

	double J = JL(dt, pL, vL, aL, AL, xL);

	double d;
	double snL = nominal_wall_thickness(xL,d);
	double AnL = nominal_area(xL,d);

	double W1L = .5*vL + 2.*pow(beta*snL/(2.*rho*AnL),.5)*(pow(AL,.25)-pow(AnL,.25));

	vnew[nx-1] = v_in;
	Anew[nx-1] = pow(.5*pow(2.*rho*Ane/(beta*sne),.5)*(-dt_act*J + W1L - .5*vnew[nx-1]) + pow(Ane,.25),4.);    
	q_in = Anew[nx-1]*vnew[nx-1];
	double dp_dA, dp_dx;
	pnew[nx-1] = pressure(x[nx-1],Anew[nx-1],dp_dA,dp_dx);
	anew[nx-1] = wave_speed(x[nx-1],Anew[nx-1]);

	double p_in = -Re*q_in + pnew[nx-1];

	return p_in;
}

//--------------------------------------------------------------
/*vector<vector<double> > moc_edge::backward_solver(vector<double> t_d, vector<double> p_d, vector<double> vfr_d)
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
}*/

//--------------------------------------------------------------
/*void moc_edge::initialization_back(vector<double> t_in, vector<double> p_in, vector<double> vfr_in)
{
	// clearing time variables
	pressure_start.clear();
   pressure_end.clear();
   velocity_start.clear();
   velocity_end.clear();
   wave_speed_start.clear();
   wave_speed_end.clear();
   total_deformation_start.clear();
   total_deformation_end.clear();
   damper_deformation_start.clear();
   damper_deformation_end.clear();
   area_start.clear();
   area_end.clear();
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
}*/

//--------------------------------------------------------------
/*double moc_edge::new_spacestep_back()
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
}*/

//--------------------------------------------------------------
/*void moc_edge::solve_back()
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
}*/

//--------------------------------------------------------------
/*void moc_edge::interpolate_back(double dx_real)
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
}*/

//--------------------------------------------------------------
/*void moc_edge::reduce_field_vectors()
{
	a.pop_back();
	d.pop_back();
	A.pop_back();
}*/

//--------------------------------------------------------------
/*double moc_edge::JA(double dt, double p, double v, double a, double d, double xp)
{
	return JR(dt, p, v, a, d, xp);
}

//--------------------------------------------------------------
double moc_edge::JB(double dt, double p, double v, double a, double d, double xp)
{
	return JL(dt, p, v, a, d, xp);
}*/

/*//--------------------------------------------------------------
void moc_edge::interpolate_hds()
{
	// interpolating along characteristic lines to dt_act
	vector<double> vq2(2*nx-2),pq2(2*nx-2),xq2(2*nx-2);
	xq2[0] = xp[0];
	vq2[0] = vp[0];
	pq2[0] = pp[0];
	xq2[2*nx-3] = xp[nx-1];
	vq2[2*nx-3] = vp[nx-1];
	pq2[2*nx-3] = pp[nx-1];
	dt[0]    = dt_act;
	dt[nx-1] = dt_act;
	for(int i=1; i<nx-1; i++)
	{
		double m = dt_act/dt[i];
		xq2[2*i-1] = m*(xp[i]-x[i-1]) + x[i-1];
		xq2[2*i]   = m*(xp[i]-x[i+1]) + x[i+1];
		vq2[2*i-1] = m*(vp[i]-v[i-1]) + v[i-1];
		vq2[2*i]   = m*(vp[i]-v[i+1]) + v[i+1];
		pq2[2*i-1] = m*(pp[i]-p[i-1]) + p[i-1];
		pq2[2*i]   = m*(pp[i]-p[i+1]) + p[i+1];
	}

	// interpolating back equidistant grid
	v[0]    = vq2[0];
	p[0]    = pq2[0];
	v[nx-1] = vq2[2*nx-3];
	p[nx-1] = pq2[2*nx-3];
	for(int i=1; i<nx-1; i++)
	{
		for(int j=0; j<2*nx-2; j++)
		{
			if(x[i] > xq2[j] && x[i] <= xq2[j+1])
			{
				v[i] = vq2[j] + (x[i]-xq2[j])/(xq2[j+1]-xq2[j])*(vq2[j+1]-vq2[j]);
				p[i] = pq2[j] + (x[i]-xq2[j])/(xq2[j+1]-xq2[j])*(pq2[j+1]-pq2[j]);
			}
		}
	}
}*/

//--------------------------------------------------------------
/*vector<double> moc_edge::boundary_start_coefficients(double dt)
{
	double xR = boundary_start_position(dt);

	// interpolating fied variables to xR
	double aR = (a[1]-a[0]) / (x[1]-x[0]) * xR + a[0];
	double vR = (v[1]-v[0]) / (x[1]-x[0]) * xR + v[0];
	double pR = (p[1]-p[0]) / (x[1]-x[0]) * xR + p[0];
	double dR = (d[1]-d[0]) / (x[1]-x[0]) * xR + d[0];

	// source term
	double J = JR(dt, pR, vR, aR, dR, xR);

	// C- characteristics constants
	double C1m = 1. / (rho*aR + Rs*A[0]) * A[0];
	double C2m = aR/(aR + Rs/rho*A[0]) * (-dt*J + vR - pR/(rho*aR)) * A[0];

	vector<double> out{C1m,C2m};
	return out;
}*/

/*//--------------------------------------------------------------
vector<double> moc_edge::boundary_end_coefficients(double dt)
{
	double xL = boundary_end_position(dt);

	// interpolating fied variables to xR
	double aL = (a[nx-2]-a[nx-1]) / (x[nx-2]-x[nx-1]) * (xL - x[nx-1]) + a[nx-1];
	double vL = (v[nx-2]-v[nx-1]) / (x[nx-2]-x[nx-1]) * (xL - x[nx-1]) + v[nx-1];
	double pL = (p[nx-2]-p[nx-1]) / (x[nx-2]-x[nx-1]) * (xL - x[nx-1]) + p[nx-1];
	double dL = (d[nx-2]-d[nx-1]) / (x[nx-2]-x[nx-1]) * (xL - x[nx-1]) + d[nx-1];

	// source term
	double J = JL(dt, pL, vL, aL, dL, xL);

	// C+ characteristics
	double C1p = -1./(rho*aL + Re*A[nx-1]) * A[nx-1];
	double C2p = aL/(aL + Re/rho*A[nx-1])*(-dt*J + vL + pL/(rho*aL)) * A[nx-1];

	vector<double> out{C1p,C2p};
	return out;
}*/


//--------------------------------------------------------------
/*vector<double> moc_edge::boundary_master_start(double dt_master)
{	
	// location of R point
	double xR = boundary_start_position(dt_master);

	// interpolating fied variables to xR
	double aR = (a[1]-a[0]) / x[1] * xR + a[0];
	double vR = (v[1]-v[0]) / x[1] * xR + v[0];
	double pR = (p[1]-p[0]) / x[1] * xR + p[0];
	double dR = (d[1]-d[0]) / x[1] * xR + d[0];

	// source term
	double J = JR(dt_master, pR, vR, aR, dR, xR);

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

	// source term
	double J = JL(dt_master, pL, vL, aL, dL, xL);

	// coefficients aq*qp + ap*pp = b
	double aq = 1./A[nx-1];
	double ap = 1./(rho*aL);
	double b  = vL + ap*pL - dt_master*J;

	vector<double> out{aq,ap,b};
	return out;
}*/

/*//--------------------------------------------------------------
double moc_edge::boundary_periferia_start(double dt, double p_out)
{
	int k=0;
	double p_s = p[0], p_old = 0., v_s = 0., delta = 0.0001;
	do
	{	
		p_old = p_s;
		double f2 = f_perif_start(p_s*(1.+delta), p_out, dt, v_s);
		double f = f_perif_start(p_s, p_out, dt, v_s);
		double fd = (f2-f)/(p_s*delta);
		p_s = p_s - f/fd;
		k++;
	}
	while(abs(p_s-p_old) > 1e-5 && k<100);

	if(k>=100)
	{
		cout << "Iteration did NOT converge at boundary_periferia_start, ID: " << ID << "\nExiting..." << endl;
		exit(-1);
	}

	// update v_s with final p_s
	f_perif_start(p_s, p_out, dt, v_s);

	// update field variables
	double ex = exp(-E2/eta2*dt_act);
	update_ith_variables(0, ex, p_s, epsz2[0], epsz[0]);

	pp[0] = p_s;
	vp[0] = v_s;
	return v_s*A[0];
}

//--------------------------------------------------------------
double moc_edge::f_perif_start(double pp, double p_out, double dt, double &v_s)
{
	double xR = boundary_start_position(dt);

	// interpolating fied variables to xR
	double aR = (a[1]-a[0]) / (x[1]-x[0]) * xR + a[0];
	double vR = (v[1]-v[0]) / (x[1]-x[0]) * xR + v[0];
	double pR = (p[1]-p[0]) / (x[1]-x[0]) * xR + p[0];
	double dR = (d[1]-d[0]) / (x[1]-x[0]) * xR + d[0];
	double AR = dR*dR*pi/4.;
	double hR = (h[1]-h[0]) / (x[1]-x[0]) * xR + h[0];
	double epszR = (epsz[1]-epsz[0]) / (x[1]-x[0]) * xR + epsz[0];
	double epsz2R = (epsz2[1]-epsz2[0]) / (x[1]-x[0]) * xR + epsz2[0];

	// source term
	double J = JR(dt, pR, vR, aR, epszR, epsz2R, dR, hR, xR, x[0]);

	v_s = vR + (pp-pR)/(rho*aR) - dt*J;

	double ex = exp(-E2/eta2*dt_act);
	double C1 = (2.*epsz[0]+1.) *dns / (E1*2.*sns*pow(epsz[0]+1.,beta)) + (2.*epsz[0]+1.)*(1.-ex) * dns / (E2*2.*sns);
	double C3 = epsz2[0]*ex;

	double f = 4.*(pp-p_out) + Rs*v_s*(dns*dns*pi*pow(C1*(pp-p0)+C3+1.,2));
	return f;
}*/

/*//--------------------------------------------------------------
double moc_edge::boundary_periferia_end(double dt, double p_out)
{
	int k=0;
	double p_e = p[nx-1], p_old = 0., v_e = 0., delta = 0.0001;
	do
	{	
		p_old = p_e;
		double f2 = f_perif_end(p_e*(1.+delta), p_out, dt, v_e);
		double f = f_perif_end(p_e, p_out, dt, v_e);
		double fd = (f2-f)/(p_e*delta);
		p_e = p_e - f/fd;
		k++;
	}
	while(abs(p_e-p_old) > 1e-5 && k<100);

	if(k>=100)
	{
		cout << "Iteration did NOT converge at boundary_periferia_end, ID: " << ID << "\nExiting..." << endl;
		exit(-1);
	}

	// update v_e with final p_e
	f_perif_end(p_e, p_out, dt, v_e);

	// update field variables
	double ex = exp(-E2/eta2*dt_act);
	update_ith_variables(nx-1, ex, p_e, epsz2[nx-1], epsz[nx-1]);

	pp[nx-1] = p_e;
	vp[nx-1] = v_e;
	return v_e*A[nx-1];
}

//--------------------------------------------------------------
double moc_edge::f_perif_end(double pp, double p_out, double dt, double &v_e)
{
	double xL = boundary_end_position(dt);

	double aL = (a[nx-2]-a[nx-1]) / (x[nx-2]-x[nx-1]) * (xL - x[nx-1]) + a[nx-1];
	double vL = (v[nx-2]-v[nx-1]) / (x[nx-2]-x[nx-1]) * (xL - x[nx-1]) + v[nx-1];
	double pL = (p[nx-2]-p[nx-1]) / (x[nx-2]-x[nx-1]) * (xL - x[nx-1]) + p[nx-1];
	double dL = (d[nx-2]-d[nx-1]) / (x[nx-2]-x[nx-1]) * (xL - x[nx-1]) + d[nx-1];
	double AL = dL*dL*pi/4.;
	double hL = (h[nx-2]-h[nx-1]) / (x[nx-2]-x[nx-1]) * (xL - x[nx-1]) + h[nx-1];
	double epszL = (epsz[nx-2]-epsz[nx-1]) / (x[nx-2]-x[nx-1]) * (xL - x[nx-1]) + epsz[nx-1];
	double epsz2L = (epsz2[nx-2]-epsz2[nx-1]) / (x[nx-2]-x[nx-1]) * (xL - x[nx-1]) + epsz2[nx-1];

	double J = JL(dt, pL, vL, aL, epszL, epsz2L, dL, hL, xL, x[nx-1]);

	v_e = vL - (pp-pL)/(rho*aL) - dt*J;

	double ex = exp(-E2/eta2*dt_act);
	double C1 = (2.*epsz[nx-1]+1.) *dne / (E1*2.*sne*pow(epsz[nx-1]+1.,beta)) + (2.*epsz[nx-1]+1.)*(1.-ex) * dne / (E2*2.*sne);
	double C3 = epsz2[nx-1]*ex;

	double f = 4.*(pp-p_out) - Re*v_e*(dne*dne*pi*pow(C1*(pp-p0)+C3+1.,2));
	return f;
}*/

/*//--------------------------------------------------------------
void moc_edge::update_variables()
{
	time.push_back(time.back() + dt_act);
	for(unsigned int i=1; i<p.size()-1; i++)
	{	
		update_ith_variables(i,p[i]);
	}
}

//--------------------------------------------------------------
void moc_edge::update_ith_variables(int i, double p_new)
{
	double dd_dp, dd_dx;
	d[i] = diameter(p_new,x[i],dd_dp,dd_dx);
	A[i] = d[i]*d[i]*pi/4.;
	a[i] = pow(d[i]/(2.*rho*dd_dp),.5);
}*/

/*//--------------------------------------------------------------
void moc_edge::interpolate()
{
	// setting back dt-s at boundaries
	dt[0]    = dt_act;
	dt[nx-1] = dt_act;

	// interpolation in space, position: xp, old time level
	vq[0]    = v[0];
	pq[0]    = p[0];
	aq[0]    = a[0];
	dq[0]    = d[0];
	vq[nx-1] = v[nx-1];
	pq[nx-1] = p[nx-1];
	aq[nx-1] = a[nx-1];
	dq[nx-1] = d[nx-1];
	cout << "q1" << endl;
	for(unsigned int i=1; i<nx-1; i++)
	{
		if(xp[i]<x[i])
		{
			double f = (xp[i]-x[i-1])/dx;
			vq[i] = v[i-1] + f*(v[i]-v[i-1]);
			pq[i] = p[i-1] + f*(p[i]-p[i-1]);
			aq[i] = a[i-1] + f*(a[i]-a[i-1]);
			dq[i] = d[i-1] + f*(d[i]-d[i-1]);
		}
		else
		{
			double f = (xp[i]-x[i])/dx;
			vq[i] = v[i] + f*(v[i+1]-v[i]);
			pq[i] = p[i] + f*(p[i+1]-p[i]);
			aq[i] = a[i] + f*(a[i+1]-a[i]);
			dq[i] = d[i] + f*(d[i+1]-d[i]);
		}
	}
	for(int i=0; i<nx; i++)
	{
		cout << "i: " << i << " p: " << pq[i] << endl;
	}
	cout << "q1" << endl;

	cout << "qt" << endl;
	// interpolation in time to new time level
	for(unsigned int i=0; i<nx; i++)
	{
		double f = dt_act/dt[i];
		vq[i] = vq[i] + f*(vp[i]-vq[i]);
		pq[i] = pq[i] + f*(pp[i]-pq[i]);
		aq[i] = aq[i] + f*(ap[i]-aq[i]);
		dq[i] = dq[i] + f*(dp[i]-dq[i]);
	}
	cout << "dt_act: " << dt_act << endl;
	for(int i=0; i<nx; i++)
	{
		cout << "i: " << i << " p: " << pq[i] << " dt: " << dt[i] << endl;
	}

	cout << "qt" << endl;

	// interpolation in space to equidistance grid
	v[0]    = vq[0];
	p[0]    = pq[0];
	a[0]    = aq[0];
	d[0]    = dq[0];
	v[nx-1] = vq[nx-1];
	p[nx-1] = pq[nx-1];
	a[nx-1] = aq[nx-1];
	d[nx-1] = dq[nx-1];
	for(unsigned int i=1; i<nx-1; i++)
	{
		if(x[i]<xp[i])
		{
			double f = (x[i]-xp[i-1])/(xp[i]-xp[i-1]);
			v[i] = vq[i-1] + f*(vq[i]-vq[i-1]);
			p[i] = pq[i-1] + f*(pq[i]-pq[i-1]);
			a[i] = aq[i-1] + f*(aq[i]-aq[i-1]);
			d[i] = dq[i-1] + f*(dq[i]-dq[i-1]);
		}
		else
		{
			double f = (x[i]-xp[i])/(xp[i+1]-xp[i]);
			v[i] = vq[i] + f*(vq[i+1]-vq[i]);
			p[i] = pq[i] + f*(pq[i+1]-pq[i]);
			a[i] = aq[i] + f*(aq[i+1]-aq[i]);
			d[i] = dq[i] + f*(dq[i+1]-dq[i]);
		}
	}

	cout << "q2" << endl;
	for(int i=0; i<nx; i++)
	{
		cout << "i: " << i << " p:" << p[i] << endl;
	}
	cout << "q2" << endl;

	// saving new time as well
	time.push_back(time.back() + dt_act);
}*/