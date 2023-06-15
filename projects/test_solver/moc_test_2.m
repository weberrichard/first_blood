
clear;

par.ds = 0.011;
par.de = 0.01;
par.ss = 0.0011;
par.se = 0.001;
par.l = 0.2;
par.E = 5e5;
par.nu = 3e-6;
par.rho = 1000;
par.nu_p = 0.5;
par.p0 = 1e5;
par.v1 = 1;

nx = 11;

dt = zeros(nx,1);
x = linspace(0,par.l,nx);
dx = par.l/(nx-1);

% initial conditions
p = par.p0*ones(nx,1);
v = zeros(nx,1);
d = zeros(nx,1);
a = zeros(nx,1);
q  = zeros(nx,1);
for i=1:nx
    [d(i),dd_dp,dd_dx] = diameter(p(i),x(i),par);
    a(i) = sqrt(1/par.rho*d(i)/2/dd_dp);
end

pnew = zeros(nx,1);
vnew = zeros(nx,1);
dnew = zeros(nx,1);
anew = zeros(nx,1);

t=[];
v1=[];
vend=[];
p1=[];
pend=[];
d1=[];
dend=[];

t_act=0;
t_end=30;
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
    
    % inner points
    for i=2:nx-1
        xL = pos_l(dt_act,v(i-1),v(i),a(i-1),a(i),dx,x(i));
        pL = p(i-1) + (xL-x(i-1))/dx*(p(i)-p(i-1));
        vL = v(i-1) + (xL-x(i-1))/dx*(v(i)-v(i-1));
        dL = d(i-1) + (xL-x(i-1))/dx*(d(i)-d(i-1));
        aL = a(i-1) + (xL-x(i-1))/dx*(a(i)-a(i-1));
        
        J_l = JL(pL,vL,aL,dL,xL,par);
        
        xR = pos_r(dt_act,v(i),v(i+1),a(i),a(i+1),dx,x(i));
        pR = p(i) + (xR-x(i))/dx*(p(i+1)-p(i));
        vR = v(i) + (xR-x(i))/dx*(v(i+1)-v(i));
        dR = d(i) + (xR-x(i))/dx*(d(i+1)-d(i));
        aR = a(i) + (xR-x(i))/dx*(a(i+1)-a(i));
        
        J_r = JR(pR,vR,aR,dR,xR,par);
        
        alpha_l = pL + par.rho*aL*vL;
        beta_r  = pR - par.rho*aR*vR;
       
        pnew(i) = .5*(alpha_l+beta_r-dt_act*(J_l+J_r));
        [dnew(i),dd_dp,dd_dx] = diameter(pnew(i),x(i),par);
        anew(i) = sqrt(1/par.rho*dnew(i)/2/dd_dp);
        vnew(i) = (alpha_l-beta_r-dt_act*(J_l-J_r))/2/par.rho/anew(i);
       
%        fprintf("P: %3i, vfr: %6.3e, a: %6.3f, d: %8.5f, p: %6.3f, v: %8.5f\n",i,vnew(i)*dnew(i)*dnew(i)*pi/4,anew(i),dnew(i),pnew(i),vnew(i));
    end
    
    % left boundary
    xR = pos_r(dt_act,v(1),v(2),a(1),a(2),dx,x(1));
    
    pR = p(1) + xR/dx*(p(2)-p(1));
    vR = v(1) + xR/dx*(v(2)-v(1));
    dR = d(1) + xR/dx*(d(2)-d(1));
    aR = a(1) + xR/dx*(a(2)-a(1));
    
    J_r = JR(pR,vR,aR,dR,xR,par);
    
    vnew(1) = par.v1;
    pnew(1) = pR + par.rho*aR*(vnew(1)-vR) - dt_act*J_r;
    [dnew(1),dd_dp,~] = diameter(pnew(1),x(1),par);
    anew(1) = sqrt(1/par.rho*dnew(1)/2/dd_dp); 
    
    % right boundary    
    xL = pos_l(dt_act,v(nx-1),v(nx),a(nx-1),a(nx),dx,x(nx));
    
    pL = p(nx-1) + (xL-x(nx-1))/dx*(p(nx)-p(nx-1));
    vL = v(nx-1) + (xL-x(nx-1))/dx*(v(nx)-v(nx-1));
    dL = d(nx-1) + (xL-x(nx-1))/dx*(d(nx)-d(nx-1));
    aL = a(nx-1) + (xL-x(nx-1))/dx*(a(nx)-a(nx-1));
    
    J_l = JL(pL,vL,aL,dL,xL,par);
    
    pnew(nx) = par.p0;
    [dnew(nx),dd_dp,~] = diameter(pnew(nx),x(nx),par);
    anew(nx) = sqrt(1/par.rho*dnew(nx)/2/dd_dp);
    vnew(nx) = vL + (pL - pnew(nx) - dt_act*J_l)/par.rho/aL;
    
    % rewriting old time level
    p = pnew;
    v = vnew;
    d = dnew;
    a = anew;
    % flow rate
    q = v.*d.^2*pi/4;
    
    % saving in time
    p1   = [p1,p(1)];
    pend = [pend,p(nx)];
    v1 = [v1,v(1)];
    vend = [vend,v(nx)];
    d1 = [d1,d(1)];
    dend = [dend,d(nx)];
    
    % plot
%     subplot(4,1,1);
%     plot(x,p);
%     subplot(4,1,2);
%     plot(x,v);
%     subplot(4,1,3);
%     plot(x,a);
%     subplot(4,1,4);
%     plot(x,d);
    asd=1;
    
end

Anx = zeros(nx,1);
A = zeros(nx,1);
for i=1:nx
   A = d(i)^2*pi/4;
   [dn,~] = nominal_diameter(x(i),par);
   Anx(i) = dn*dn*pi/4;
end

subplot(4,1,1);
plot(x,p);
subplot(4,1,2);
plot(x,v);
hold on;
plot(x,v+a,'r');
plot(x,v-a,'r');
subplot(4,1,3);
plot(x,A);
hold on;
plot(x,Anx,'k');
subplot(4,1,4);
plot(x,q);

function JL = JL(p,v,a,d,x,par)
    [~,dd_dp,dd_dx] = diameter(p,x,par);
    JL = par.rho*a*32*par.nu/d/d*v + v*dd_dx/dd_dp;
%     JL = 0.02/2/d*v* abs(v)*par.rho*a + v*dd_dx/dd_dp;
end

function JR = JR(p,v,a,d,x,par)
    [~,dd_dp,dd_dx] = diameter(p,x,par);
    JR = -par.rho*a*32*par.nu/d/d*v + v*dd_dx/dd_dp;
%     JR = -0.02/2/d*v* abs(v)*par.rho*a + v*dd_dx/dd_dp;
end

function [d,dd_dp,dd_dx] = diameter(p,x,par)
    [dn,dn_dx] = nominal_diameter(x,par);
    [sn,sn_dx] = nominal_wall(x,par);
    
    d     = dn + (1-par.nu_p^2)*dn^2/2/sn/par.E*(p-par.p0);
    dd_dp = (1-par.nu_p^2)*dn^2/2/sn/par.E;
    dd_dx = dn_dx + (p-par.p0)*(1-par.nu_p^2)/2/par.E*(2*dn/sn*dn_dx-dn^2/sn^2*sn_dx);
end

function [dn,dn_dx] = nominal_diameter(x,par)
    dn = par.ds + x/par.l*(par.de-par.ds);
    dn_dx = (par.de-par.ds)/par.l;
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

