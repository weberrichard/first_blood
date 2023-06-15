
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
par.CFL = 0.9;

nx = 11;

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
Aref1 = max(Anx);
Aref2 = min(Anx);
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

U=zeros(2,nx);
F=zeros(2,nx);
S=zeros(2,nx);
Us=zeros(2,nx);
Fs=zeros(2,nx);
Ss=zeros(2,nx);

t_act=0;
t_end=4.9;
while(t_act<t_end)
   
    % finding dt-s
    dt(1) = dx/(a(2)-v(2));
    for i=2:nx-1
        dt(i) = min(dx/(a(i-1)+v(i-1)),dx/(a(i+1)-v(i+1)));
    end
    dt(nx) = dx/(a(nx-1)+v(nx-1));
    dt_act = par.CFL*min(dt);
    t_act = t_act + dt_act;
    t = [t,t_act];
    
    % inner points
    Sp = zeros(nx,1);
    for i=1:nx
        [~,~,dp_dx] = cross_section(Anew(i),x(i),par);
        Sp(i) = dp_dx;
    end
    U = [A';v'];
    F = [(A.*v)';(v.^2/2+p/par.rho)'];
    S = [zeros(1,nx);(-8*pi*par.nu./A.*v)'];

    Us(:,1) = U(:,1) - dt_act/dx*(F(:,2)-F(:,1)) + dt_act*S(:,1);
    As = Us(1,1);
    vs = Us(2,1);
    [ps,~,dp_dxs] = cross_section(As,x(1),par);
    Fs(:,1) = [As*vs;vs^2/2+ps/par.rho];
    for i=2:nx-1
        Us(:,i) = U(:,i) - dt_act/dx*(F(:,i+1)-F(:,i)) + dt_act*S(:,i);

        As = Us(1,i);
        vs = Us(2,i);
        [ps,~,dp_dxs] = cross_section(As,x(i),par);
        Fs(:,i) = [As*vs;vs^2/2+ps/par.rho];
        Ss(:,i) = [0;-8*pi*par.nu/As*vs];

        U(:,i) = .5*(U(:,i)+Us(:,i)) - dt_act/2/dx*(Fs(:,i)-Fs(:,i-1)) + dt_act/2*Ss(:,i);
        Anew(i) = U(1,i);
        vnew(i) = U(2,i);
        [pnew(i),dp_dA,~] = cross_section(Anew(i),x(i),par);
        anew(i) = sqrt(Anew(i)/par.rho*dp_dA);
    end

    % left boundary - inlet velocity
    xR = pos_r(dt_act,v(1),v(2),a(1),a(2),dx,x(1));

    pR = p(1) + xR/dx*(p(2)-p(1));
    vR = v(1) + xR/dx*(v(2)-v(1));
    AR = A(1) + xR/dx*(A(2)-A(1));

    J_r = JR(AR,vR,xR,par);

    [dnR,~] = nominal_diameter(xR,par);
    [snR,~] = nominal_wall(xR,par);
    AnR = dnR*dnR*pi/4;
    
    W2R = -4* sqrt(par.beta*snR/(2*par.rho*AnR)) *(AR^.25-Aref1^.25) + vR;

    [dn,~] = nominal_diameter(x(1),par);
    [sn,~] = nominal_wall(x(1),par);
    An = dn*dn*pi/4;

    vnew(1) = par.v1;
    Anew(1) = 1/64*(par.rho*An/(par.beta*sn))^2 * (vnew(1) - W2R + dt_act*J_r + 4* sqrt(par.beta*sn/(2*par.rho*An))*Aref1^.25)^4;
    [pnew(1),dp_dA,dp_dx] = cross_section(Anew(1),x(1),par);
    anew(1) = sqrt(Anew(1)/par.rho*dp_dA);

    % right boundary - outlet pressure
    xL = pos_l(dt_act,v(nx-1),v(nx),a(nx-1),a(nx),dx,x(nx));

    pL = p(nx-1) + (xL-x(nx-1))/dx*(p(nx)-p(nx-1));
    vL = v(nx-1) + (xL-x(nx-1))/dx*(v(nx)-v(nx-1));
    AL = A(nx-1) + (xL-x(nx-1))/dx*(A(nx)-A(nx-1));

    J_l = JL(AL,vL,xL,par);

    [dnL,~] = nominal_diameter(xL,par);
    [snL,~] = nominal_wall(xL,par);
    AnL = dnL*dnL*pi/4;

    W1L =  4* sqrt(par.beta*snL/(2*par.rho*AnL)) *(AL^.25-Aref2^.25) + vL;

    [dn,~] = nominal_diameter(x(nx),par);
    [sn,~] = nominal_wall(x(nx),par);
    An = dn*dn*pi/4;

    pnew(nx) = par.p0;
    Anew(nx) = ( (pnew(nx)-par.p0) *An/(par.beta*sn) + sqrt(An) )^2;
    vnew(nx) = W1L - dt_act*J_l + 4* sqrt(par.beta*sn/(2*par.rho*An))*(Aref2^.25-Anew(nx)^.25);
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

%     for i=1:nx
%         fprintf("i: %3i, A: %6.3e, v: %6.3f, p: %8.1f, a: %6.3f\n",i,Anew(i),vnew(i),pnew(i),anew(i));
%     end
%     fprintf("\n");

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
plot(x,p-par.p0);
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

function J_l = JL(A,v,x,par)
    [~,~,dp_dx] = cross_section(A,x,par);
    J_l = 8*pi*par.nu/A*v + 1/par.rho*dp_dx;
end

function J_r = JR(A,v,x,par)
    [~,~,dp_dx] = cross_section(A,x,par);
    J_r = 8*pi*par.nu/A*v + 1/par.rho*dp_dx;
end

function [p,dp_dA,dp_dx] = cross_section(A,x,par)
    [dn,dn_dx] = nominal_diameter(x,par);
    [sn,sn_dx] = nominal_wall(x,par);
    An = dn*dn*pi/4;
    An_dx = dn*pi/2*dn_dx;
    p = par.p0 + par.beta*sn/An*(sqrt(A)-sqrt(An));
    dp_dx = par.beta/An*( (sqrt(A)-sqrt(An))*(sn_dx - sn/An*An_dx) - .5*sn/sqrt(An)*An_dx );
%     dp_dx = 0;
    dp_dA = par.beta*sn/An*.5/sqrt(A);
end

function [dn,dn_dx] = nominal_diameter(x,par)
%     dn = par.ds + x/par.l*(par.de-par.ds);
%     dn_dx = (par.de-par.ds)/par.l;
    dn = sqrt(par.ds^2+x/par.l*(par.de^2-par.ds^2));
    dn_dx = 1/2/dn*(par.de^2-par.ds^2)/par.l;
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