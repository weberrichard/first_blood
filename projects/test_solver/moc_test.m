clear;

par.ds = 0.01;
par.de = 0.01;
par.ss = 0.001;
par.se = 0.001;
par.l = 10.2;
par.E = 5e5;
par.nu = 3e-6;
par.rho = 1055;
par.nu_p = 0.5;
par.p0 = 1e5;

nx = 51;

v0 = 1;
pn = 1e5;

pp = zeros(nx,1);
vp = zeros(nx,1);
ap = zeros(nx,1);
dp = zeros(nx,1);
q  = zeros(nx,1);

dt = zeros(nx,1);
x = linspace(0,par.l,nx);
dx = par.l/(nx-1);
xp = zeros(nx,1);

% initial conditions
p = par.p0*ones(nx,1);
v = zeros(nx,1);
d = zeros(nx,1);
a = zeros(nx,1);
for i=1:nx
    %[d(i),dd_dp,dd_dx] = diameter(p(i),x(i),par);
    %a(i) = sqrt(1/par.rho*d(i)/2/dd_dp);
    d(i) = par.de;
    a(i) = 10;
end

t=[];
v1=[];
vend=[];
p1=[];
pend=[];
d1=[];
dend=[];

t_act=0;
t_end=100;
while(t_act<t_end)
   
    % finding dt-s
    dt(1) = dx/(a(2)-v(2));
    for i=2:nx-1
       dt(i) = 2*dx/(a(i-1)+v(i-1)+a(i+1)-v(i+1));
    end
    dt(nx) = dx/(a(nx-1)+v(nx-1));
    dt_act = min(dt);
    t_act = t_act + dt_act;
    t = [t,t_act];
    
    % inner points
    for i=2:nx-1
       xp(i) = x(i-1)+dt(i)*(a(i-1)+v(i-1));
       
       alpha_l = p(i-1) + par.rho*a(i-1)*v(i-1);
       beta_r  = p(i+1) - par.rho*a(i+1)*v(i+1);
       
       J_l = JL(p(i-1),v(i-1),a(i-1),d(i-1),x(i-1),par);
       J_r = JR(p(i+1),v(i+1),a(i+1),d(i+1),x(i+1),par);
       
       pp(i) = .5*(alpha_l+beta_r-dt(i)*(J_l+J_r));
%        [dp(i),dd_dp,dd_dx] = diameter(pp(i),xp(i),par);
%        ap(i) = sqrt(1/par.rho*dp(i)/2/dd_dp);
%        vp(i) = (alpha_l-beta_r-dt(i)*(J_l-J_r))/2/par.rho/ap(i);
       dp(i) = par.de;
       ap(i) = 10;
       vp(i) = (alpha_l-beta_r-dt(i)*(J_l-J_r))/2/par.rho/ap(i);
       
       fprintf("P: %3i, vfr: %6.3e, a: %6.3f, d: %8.5f, p: %6.3f, v: %8.5f\n",i,v(i)*dp(i)*dp(i)*pi/4,ap(i),dp(i),pp(i),vp(i));
    end
    
    % left boundary
    vp(1) = v0;
    
    dv = (v(2)-v(1)) / dx;
	da = (a(2)-a(1)) / dx;
	xR = dt(1)*(a(1)-v(1)) / (1.+dt(1)*dv-dt(1)*da);
    
    pR = p(1) + xR/dx*(p(2)-p(1));
    vR = v(1) + xR/dx*(v(2)-v(1));
    dR = d(1) + xR/dx*(d(2)-d(1));
    aR = a(1) + xR/dx*(a(2)-a(1));
    
    J_r = JR(pR,vR,aR,dR,xR,par);
    
    pp(1) = pR + par.rho*aR*(vp(1)-vR) - dt(1)*J_r;
%     [dp(1),dd_dp,~] = diameter(pp(1),xR,par);
%     ap(1) = sqrt(1/par.rho*dp(1)/2/dd_dp);   
    dp(1) = par.de;
    ap(1) = 10;
    
    % right boundary    
    dv = (v(nx)-v(nx-1)) / dx;
	da = (a(nx)-a(nx-1)) / dx;
	xL = x(nx) - dt(nx)*(v(nx)+a(nx)) / (1.+dt(nx)*da+dt(nx)*dv);
    
    pL = p(nx-1) + (xL-x(nx-1))/dx*(p(nx)-p(nx-1));
    vL = v(nx-1) + (xL-x(nx-1))/dx*(v(nx)-v(nx-1));
    dL = d(nx-1) + (xL-x(nx-1))/dx*(d(nx)-d(nx-1));
    aL = a(nx-1) + (xL-x(nx-1))/dx*(a(nx)-a(nx-1));
    
    J_l = JL(pL,vL,aL,dL,xL,par);
    
    pp(nx) = pn;
%     [dp(nx),dd_dp,~] = diameter(pp(nx),xL,par);
%     ap(nx) = sqrt(1/par.rho*dp(nx)/2/dd_dp);
    dp(nx) = par.de;
    ap(nx) = 10;
    vp(nx) = vL + (pL - pp(nx) - dt(nx)*J_l)/par.rho/aL;
    xp(nx) = par.l;
    
    % interpolation
    dt(1) = dt_act;
    dt(nx) = dt_act;
    pq = interp1(x,p,xp);
    vq = interp1(x,v,xp);
    aq = interp1(x,a,xp);
    dq = interp1(x,d,xp);
    
    for i=1:nx
       pq(i) = pq(i) + dt_act/dt(i)*(pp(i)-pq(i));
       vq(i) = vq(i) + dt_act/dt(i)*(vp(i)-vq(i));
       aq(i) = aq(i) + dt_act/dt(i)*(ap(i)-aq(i));
       dq(i) = dq(i) + dt_act/dt(i)*(dp(i)-dq(i));
    end
    
    p = interp1(xp,pq,x);
    v = interp1(xp,vq,x);
    a = interp1(xp,aq,x);
    d = interp1(xp,dq,x);
    
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
    
    fprintf("\n");
    for i=1:nx
        fprintf("%3i, vfr: %6.3e, a: %6.3f, d: %8.5f, p: %6.3f, v: %8.5f\n",i,q(i),a(i),d(i),p(i),v(i));
    end
    asd=1;
    
end

subplot(4,1,1);
plot(x,p);
subplot(4,1,2);
plot(x,v);
subplot(4,1,3);
plot(x,a);
subplot(4,1,4);
plot(x,d);

function JL = JL(p,v,a,d,x,par)
    %[~,dd_dp,dd_dx] = diameter(p,x,par);
    JL = par.rho*a*32*par.nu/d/d*v; % + v*dd_dx/dd_dp;
end

function JR = JR(p,v,a,d,x,par)
    %[~,dd_dp,dd_dx] = diameter(p,x,par);
    JR = -par.rho*a*32*par.nu/d/d*v; % + v*dd_dx/dd_dp;
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

% function xL = pos_l(dt,v1,v2,a1,a2,dx,x)
%     dv = (v2-v1) / dx;
% 	da = (a2-a1) / dx;
% 	xL = x - dt*(v2+a2) / (1.+dt*da+dt*dv);
% end
% 
% function xR = pos_r(dt,v1,v2,a1,a2,dx)
%     dv = (v2-v1) / dx;
% 	da = (a2-a1) / dx;
% 	xR = x + dt*(a1-v1) / (1.+dt1*dv-dt1*da);
% end

