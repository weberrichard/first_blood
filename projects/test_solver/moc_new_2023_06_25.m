
clear;

par.ds = 0.015;
par.de = 0.01;
par.As = par.ds*par.ds*pi/4;
par.Ae = par.de*par.de*pi/4;
par.ss = 0.0010;
par.se = 0.0010;
par.l = 0.1;
par.E = 5e5;
par.nu = 3e-6;
par.rho = 1055;
par.nu_p = 0.5;
par.p0 = 1e5;
par.v1 = 1;
par.beta = sqrt(pi)*par.E/(1-par.nu_p^2);

nx = 41;

dt = zeros(nx,1);
x = linspace(0,par.l,nx);
dx = par.l/(nx-1);

% initial conditions
p = par.p0*ones(nx,1);
v = zeros(nx,1);
A = zeros(nx,1);
a = zeros(nx,1);
q  = zeros(nx,1);
dnx = zeros(nx,1);
for i=1:nx
    [dn,~] = nominal_diameter(x(i),par);
    A(i) = dn*dn*pi/4;
    dnx(i) = dn;
    [p(i),dp_dA,dp_dx] = cross_section(A(i),x(i),par);
    a(i) = sqrt(A(i)/par.rho*dp_dA);
end

Vmax = sum(A(1:nx-1)*dx);
Vmin = sum(A(2:nx)*dx);
Vmean = (Vmax+Vmin)*.5;
Amean = Vmean/par.l;

Anx = A;
Aref1 = Amean;
Aref2 = Amean;
pnew = zeros(nx,1);
vnew = zeros(nx,1);
Anew = zeros(nx,1);
anew = zeros(nx,1);

t=[];
v1=[];
vend=[];
p1=[];
pend=[];
A1=[];
Aend=[];

t_act=0;
t_end=20;
while(t_act<t_end)
   
    % finding dt-s
    dt(1) = dx/(a(2)-v(2));
    for i=2:nx-1
        dt(i) = min(dx/(a(i-1)+v(i-1)),dx/(a(i+1)-v(i+1)));
    end
    dt(nx) = dx/(a(nx-1)+v(nx-1));
    dt_act = min(dt);
    t_act = t_act + dt_act;
    t = [t,t_act];

    for i=2:nx-1
        xL = pos_l(dt_act,v(i-1),v(i),a(i-1),a(i),dx,x(i));
        vL = v(i-1) + (xL-x(i-1))/dx*(v(i)-v(i-1));
        AL = A(i-1) + (xL-x(i-1))/dx*(A(i)-A(i-1));
        aL = a(i-1) + (xL-x(i-1))/dx*(a(i)-a(i-1));
        
        J_l = JL(AL,vL,xL,par,aL);

        [dnL,~] = nominal_diameter(xL,par);
        [snL,~] = nominal_wall(xL,par);
        AnL = dnL*dnL*pi/4;
        
        xR = pos_r(dt_act,v(i),v(i+1),a(i),a(i+1),dx,x(i));
        vR = v(i) + (xR-x(i))/dx*(v(i+1)-v(i));
        AR = A(i) + (xR-x(i))/dx*(A(i+1)-A(i));
        aR = a(i) + (xR-x(i))/dx*(a(i+1)-a(i));
        
        J_r = JR(AR,vR,xR,par,aR);
           
        [dnR,~] = nominal_diameter(xR,par);
        [snR,~] = nominal_wall(xR,par);
        AnR = dnR*dnR*pi/4;
        
        [dn,~] = nominal_diameter(x(i),par);
        [sn,~] = nominal_wall(x(i),par);
        An = dn*dn*pi/4;
        
        W1L =  4* sqrt(par.beta*snL/(2*par.rho*AnL)) *(AL^.25-AnL^.25) + vL;
        W2R = -4* sqrt(par.beta*snR/(2*par.rho*AnR)) *(AR^.25-AnR^.25) + vR;
        
        vnew(i) = .5*(W1L+W2R-dt_act*(J_l+J_r) + 4* sqrt(par.beta*sn/(2*par.rho*An))*(An^.25-An^.25));
        Anew(i) = 1/1024 * (par.rho*An/(par.beta*sn))^2 * (W1L-W2R-dt_act*(J_l-J_r) ...
            + 4* sqrt(par.beta*sn/(2*par.rho*An))*(An^.25+An^.25))^4;
        [pnew(i),dp_dA,~] = cross_section(Anew(i),x(i),par);
        anew(i) = sqrt(Anew(i)/par.rho*dp_dA);
       
%        fprintf("P: %3i, vfr: %6.3e, a: %6.3f, d: %8.5f, p: %6.3f, v: %8.5f\n",i,vnew(i)*dnew(i)*dnew(i)*pi/4,anew(i),dnew(i),pnew(i),vnew(i));
    end

    % left boundary - inlet velocity
    xR = pos_r(dt_act,v(1),v(2),a(1),a(2),dx,x(1));

    pR = p(1) + xR/dx*(p(2)-p(1));
    vR = v(1) + xR/dx*(v(2)-v(1));
    AR = A(1) + xR/dx*(A(2)-A(1));
    aR = a(1) + xR/dx*(a(2)-a(1));

    J_r = JR(AR,vR,xR,par,aR);

    [dnR,~] = nominal_diameter(xR,par);
    [snR,~] = nominal_wall(xR,par);
    AnR = dnR*dnR*pi/4;
    W2R = -4* sqrt(par.beta*snR/(2*par.rho*AnR)) *(AR^.25-AnR^.25) + vR;

    [dn,~] = nominal_diameter(x(1),par);
    [sn,~] = nominal_wall(x(1),par);
    An = dn*dn*pi/4;

    vnew(1) = par.v1;
    Anew(1) = 1/64*(par.rho*An/(par.beta*sn))^2 * (vnew(1) - W2R + dt_act*J_r + 4* sqrt(par.beta*sn/(2*par.rho*An))*An^.25)^4;
    [pnew(1),dp_dA,dp_dx] = cross_section(Anew(1),x(1),par);
    anew(1) = sqrt(Anew(1)/par.rho*dp_dA);

    % right boundary - outlet pressure
    xL = pos_l(dt_act,v(nx-1),v(nx),a(nx-1),a(nx),dx,x(nx));

    pL = p(nx-1) + (xL-x(nx-1))/dx*(p(nx)-p(nx-1));
    vL = v(nx-1) + (xL-x(nx-1))/dx*(v(nx)-v(nx-1));
    AL = A(nx-1) + (xL-x(nx-1))/dx*(A(nx)-A(nx-1));
    aL = a(nx-1) + (xL-x(nx-1))/dx*(a(nx)-a(nx-1));

    J_l = JL(AL,vL,xL,par,aL);

    [dnL,~] = nominal_diameter(xL,par);
    [snL,~] = nominal_wall(xL,par);
    AnL = dnL*dnL*pi/4;
    W1L =  4* sqrt(par.beta*snL/(2*par.rho*AnL)) *(AL^.25-AnL^.25) + vL;

    [dn,~] = nominal_diameter(x(nx),par);
    [sn,~] = nominal_wall(x(nx),par);
    An = dn*dn*pi/4;

    pnew(nx) = par.p0;
    Anew(nx) = ( (pnew(nx)-par.p0) *An/(par.beta*sn) + sqrt(An) )^2;
    vnew(nx) = W1L - dt_act*J_l + 4* sqrt(par.beta*sn/(2*par.rho*An))*(An^.25-Anew(nx)^.25);
    [~,dp_dA,~] = cross_section(Anew(nx),x(nx),par);
    anew(nx) = sqrt(Anew(nx)/par.rho*dp_dA);
 
    % rewriting old time level
    p = pnew;
    v = vnew;
    A = Anew;
    a = anew;
    % flow rate
    q = v.*A;
    d = sqrt(A*4/pi);
    
    % saving in time
    p1   = [p1,p(1)];
    pend = [pend,p(nx)];
    v1 = [v1,v(1)];
    vend = [vend,v(nx)];
    A1 = [A1,A(1)];
    Aend = [Aend,A(nx)];
    
    % plot
%     figure(1);
%     subplot(4,1,1);
%     plot(x,p);
%     subplot(4,1,2);
%     cla;
%     plot(x,v);
%     hold on
%     plot(x,a,'r');
%     plot(x,-a,'r');
%     subplot(4,1,3);
%     cla;
%     plot(x,A);
%     hold on;
%     plot(x,Anx,'k');
%     subplot(4,1,4);
%     cla;
%     plot(x,d);
%     hold on;
%     plot(x,dnx,'k');
    asd=1;
end

figure;
subplot(4,1,1);
plot(x,p);
ylabel('p [Pa]');
subplot(4,1,2);
plot(x,v);
hold on;
plot(x,a,'r');
plot(x,-a,'r');
ylabel('v [m/s]');
subplot(4,1,3);
plot(x,A);
hold on;
plot(x,Anx,'k');
ylabel('A [m2]');
subplot(4,1,4);
plot(x,q);
xlabel('x [m]');
ylabel('q [m3/s');

figure;
subplot(2,1,1);
plot(t,(p1-par.p0)/133.36);
ylabel('p [mmHg]');
subplot(2,1,2);
plot(t,vend);
xlabel('t [s]');
ylabel('v [m/s]');

function J_l = JL(A,v,x,par,a)
    [~,~,dp_dx] = cross_section(A,x,par);
    [dn,dn_dx] = nominal_diameter(x,par);
    [sn,sn_dx] = nominal_wall(x,par);
    An = dn*dn*pi/4;
    An_dx = dn*pi/2*dn_dx;
    C1 = -2*sqrt(par.beta/2/par.rho)*((0.5*sqrt(1/An/sn)*sn_dx-0.5*sqrt(sn/An^3)*An_dx)*(A^0.25-An^0.25)-sqrt(sn/An)*0.25*An^(-0.75)*An_dx);
    J_l = 8*pi*par.nu/A*v + 1/par.rho*dp_dx + 2*(v+a)*C1;
end

function J_r = JR(A,v,x,par,a)
    [~,~,dp_dx] = cross_section(A,x,par);
    [dn,dn_dx] = nominal_diameter(x,par);
    [sn,sn_dx] = nominal_wall(x,par);
    An = dn*dn*pi/4;
    An_dx = dn*pi/2*dn_dx;
    C2 = 2*sqrt(par.beta/2/par.rho)*((0.5*sqrt(1/An/sn)*sn_dx-0.5*sqrt(sn/An^3)*An_dx)*(A^0.25-An^0.25)-sqrt(sn/An)*0.25*An^(-0.75)*An_dx);
    J_r = 8*pi*par.nu/A*v + 1/par.rho*dp_dx + 2*(v-a)*C2;
end

function [p,dp_dA,dp_dx] = cross_section(A,x,par)
    [dn,dn_dx] = nominal_diameter(x,par);
    [sn,sn_dx] = nominal_wall(x,par);
    An = dn*dn*pi/4;
    An_dx = dn*pi/2*dn_dx;
    p = par.p0 + par.beta*sn/An*(sqrt(A)-sqrt(An));
    dp_dx = par.beta/An*( (sqrt(A)-sqrt(An))*(sn_dx - sn/An*An_dx) - .5*sn/sqrt(An)*An_dx );
       %+ 1/An*(p-par.p0+par.beta*sn/2/sqrt(An))*An_dx;
    % dp_dx = 0;
    dp_dA = par.beta*sn/An*.5/sqrt(A);
end

function [dn,dn_dx] = nominal_diameter(x,par)
    dn = par.ds + x/par.l*(par.de-par.ds);
    dn_dx = (par.de-par.ds)/par.l;
%     dn = sqrt(par.ds^2+x/par.l*(par.de^2-par.ds^2));
%     dn_dx = 1/2/dn*(par.de^2-par.ds^2)/par.l;
end

function [sn,sn_dx] = nominal_wall(x,par)
    sn = par.ss + x/par.l*(par.se-par.ss);
    sn_dx = (par.se-par.ss)/par.l;
end

function xL = pos_l(dt,v1,v2,a1,a2,dx,x)
    dv = (v1-v2) / dx;
	da = (a1-a2) / dx;
	xL = x - dt*(v2+a2) / (1.-dt*(da+dv));
end

function xR = pos_r(dt,v2,v3,a2,a3,dx,x)
    dv = (v3-v2) / dx;
	da = (a3-a2) / dx;
	xR = x + dt*(a2-v2) / (1.-dt*(da-dv));
end