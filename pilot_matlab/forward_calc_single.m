clear;

% gemoetrical parameters
L  = 0.109193669932178;
Dn = 0.026;
dn = 0.0026;

% material parameters
E1    = 751821.278668506;
E2    =  769729.996566565;
eta2  =   13000;
beta  = 2.0;

% physical parameters
g    = 9.81;
nu   = 3e-6;
rho  = 1050;
p0   = 1e5;
mmHg_Pa = 133.322368;

% boundary conditions
A = readmatrix('heart_cycle.txt');
t_in_0 = A(:,1);
t_in = [];
for i=1:10
    if(i==1)
        t_in = t_in_0;
    else
        t_in = [t_in;0.001+t_in(end)+t_in_0];
    end
end
p_in = mmHg_Pa*A(:,2) + p0;
p_in = repmat(p_in,[10,1]);
% p_out = mmHg_Pa*(A(:,2)-10) + p0;
p_out = p0;
dzeta_out = 100000; 

% numerical pars
CFL = 1;
Nx  = 29;
T  = 0.00373; % wim in sec

x = linspace(0,L,Nx);
dx = x(2)-x(1);
t = 0;

% initial conditions
h = zeros(Nx,1);
v = zeros(Nx,1);
p = p0*ones(Nx,1);
% p = linspace(p_in(1),p_out(1),Nx);
epsz2 = zeros(Nx,1);
epsz = zeros(Nx,1);
D = Dn*ones(Nx,1);
a = sqrt(E1*dn/Dn/rho)*ones(Nx,1);

% temp variables
vp = zeros(Nx,1);
pp = zeros(Nx,1);
dt = zeros(Nx,1);
xp = zeros(Nx,1);
vq = zeros(Nx,1);
pq = zeros(Nx,1);

% saving time vars
Nt = round(1.2*T*10/dx);
pt = zeros(Nx,Nt);
vt = zeros(Nx,Nt);
at = zeros(Nx,Nt);
Dt = zeros(Nx,Nt);
epszt = zeros(Nx,Nt);
epsz2t = zeros(Nx,Nt);
pt(:,1) = p;
vt(:,1) = v;
at(:,1) = a;
Dt(:,1) = D;
epszt(:,1) = epsz;
epsz2t(:,1) = epsz2;

figure('position',[400,200,1000,750]);

j = 1;
while(t(j)<T)
    for i=2:Nx-1
       % calculating the time and location of P point
       dt(i) =  2*dx / (v(i-1)+a(i-1)-v(i+1)+a(i+1));
       xp(i) = x(i-1)+dt(i)*(v(i-1)+a(i-1));
       
       AL = 1/eta2 * (p(i-1)-p0) * Dn / (2*dn)*(2*epsz(i-1)+1) * exp(-E2/eta2*dt(i)) ...
           - E2/eta2*epsz2(i-1)*exp(-E2/eta2*dt(i));
       AR = 1/eta2 * (p(i+1)-p0) * Dn / (2*dn)*(2*epsz(i+1)+1) * exp(-E2/eta2*dt(i)) ...
           - E2/eta2*epsz2(i+1)*exp(-E2/eta2*dt(i));
       
       
       hP = (xp(i)-x(i-1))/(2*dx)*h(i-1) + (xp(i)-x(i+1))/(-2*dx)*h(i+1);
       JL = g*(hP-h(i-1))/(xp(i)-x(i-1)) + 32*nu/D(i-1)^2*v(i-1) + 2*a(i-1)/(2*epsz(i-1)+1)*AL;
       JR = g*(hP-h(i+1))/(xp(i)-x(i+1)) + 32*nu/D(i+1)^2*v(i+1) - 2*a(i+1)/(2*epsz(i+1)+1)*AR;
       
       % new velocity and pressure at P point
       vp(i) = (p(i-1)-p(i+1) + rho*a(i-1)*v(i-1) + rho*a(i+1)*v(i+1) - dt(i)*rho*a(i-1)*JL - dt(i)*rho*a(i+1)*JR) ...
           / (rho*a(i-1)+rho*a(i+1));
       pp(i) = p(i+1) + (vp(i)-v(i+1))*rho*a(i+1) + rho*a(i+1)*dt(i)*JR;
    end

    
    % left boundary conditions
    dt(1) = abs(dx/(v(2)-a(2)));
    xp(1) = 0;
    pp(1) = interp1(t_in,p_in,t(j)+dt(1));
    AR = 1/eta2 * (p(2)-p0) * Dn / (2*dn)*(2*epsz2(2)+1) * exp(-E2/eta2*dt(1)) ...
        - E2/eta2*epsz2(2)*exp(-E2/eta2*dt(1));
    
    hP = h(1);
    JR = g*(hP-h(2))/dx + 32*nu/D(2)^2*v(2) - 2*a(2)/(2*epsz(2)+1)*AR;
    vp(1) = v(2)+1/rho/a(2)*(pp(1)-p(2))-dt(1)*JR;  
    
    % right boundary conditions
    dt(Nx) = dx/(v(Nx-1)+a(Nx-1));
    xp(Nx) = L;
    A_out = D(Nx)*D(Nx)*pi/4;

    AL = 1/eta2 * (p(Nx-1)-p0) * Dn / (2*dn)*(2*epsz2(Nx-1)+1) * exp(-E2/eta2*dt(Nx)) ...
        - E2/eta2*epsz2(Nx-1)*exp(-E2/eta2*dt(Nx));
    
    hP = h(Nx);
    JL =  g*(hP-h(Nx-1))/dx + 32*nu/D(Nx-1)^2*v(Nx-1) + 2*a(Nx-1)/(2*epsz(Nx-1)+1)*AL;
    vp(Nx) = (v(Nx-1) - dt(Nx)*JL + 1/rho/a(Nx-1)*(p(Nx-1)-p0)) / (1 + 1/a(Nx-1)*dzeta_out*A_out);
    pp(Nx) = p0 + dzeta_out*rho*A_out*vp(Nx);
    
    % calculating the real time step
    dt_real = CFL*min(dt);
    t = [t,t(end)+dt_real];
    
    % interpolation in space to Q point %
    % second order
    vq(1) = v(1); pq(1) = p(1);
    vq(Nx) = v(Nx); pq(Nx) = p(Nx);
    for i=2:Nx-1
        a1 =  (xp(i)-x(i))*(xp(i)-x(i+1)) / (2*dx*dx);
        a2 =  (xp(i)-x(i-1))*(xp(i)-x(i+1)) / (-dx*dx);
        a3 =  (xp(i)-x(i-1))*(xp(i)-x(i)) / (2*dx*dx);

        vq(i) = (a1*v(i-1)+a2*v(i)+a3*v(i+1));
        pq(i) = (a1*p(i-1)+a2*p(i)+a3*p(i+1));
    end
    % first order
    vq(1) = v(1); pq(1) = p(1);
    vq(Nx) = v(Nx); pq(Nx) = p(Nx);
    for i=2:Nx-1
       if(xp(i)<x(i))
           vq(i) = v(i-1) + (xp(i)-x(i-1))/dx*(v(i)-v(i-1));
           pq(i) = p(i-1) + (xp(i)-x(i-1))/dx*(p(i)-p(i-1));
       else
           vq(i) = v(i) + (xp(i)-x(i))/dx*(v(i+1)-v(i));
           pq(i) = p(i) + (xp(i)-x(i))/dx*(p(i+1)-p(i));
       end
    end
    
    % first order with matlab
%     vq = interp1(x,v,xp);
%     pq = interp1(x,p,xp);
    % in time % 
    for i=1:Nx 
       vq(i) = vq(i) + dt_real/dt(i)*(vp(i)-vq(i));
       pq(i) = pq(i) + dt_real/dt(i)*(pp(i)-pq(i)); 
    end
    % in space %
    % second-order
%     v1(1) = vq(1);  p1(1) = pq(1);
%     v1(Nx) = vq(Nx); p1(Nx) = pq(Nx);
%     for i=2:Nx-1
%         a1 = (x(i)-xp(i))*(x(i)-xp(i+1)) / ((xp(i-1)-xp(i))*(xp(i-1)-xp(i+1)));
%         a2 = (x(i)-xp(i-1))*(x(i)-xp(i+1)) / ((xp(i)-xp(i-1))*(xp(i)-xp(i+1)));
%         a3 = (x(i)-xp(i-1))*(x(i)-xp(i)) / ((xp(i+1)-xp(i-1))*(xp(i+1)-xp(i)));
%         v1(i) = a1*vq(i-1) + a2*vq(i) + a3*vq(i+1);
%         p1(i) = a1*pq(i-1) + a2*pq(i) + a3*pq(i+1);
%     end
    
    % linear, manually
    v(1) = vq(1);  p(1) = pq(1);
    v(Nx) = vq(Nx); p(Nx) = pq(Nx);
    for i=2:Nx-1
        if(x(i)<xp(i))
            v(i) = vq(i-1) + (x(i)-xp(i-1))/(xp(i)-xp(i-1))*(vq(i)-vq(i-1));
            p(i) = pq(i-1) + (x(i)-xp(i-1))/(xp(i)-xp(i-1))*(pq(i)-pq(i-1));
        else
            v(i) = vq(i) + (x(i)-xp(i))/(xp(i+1)-xp(i))*(vq(i+1)-vq(i));
            p(i) = pq(i) + (x(i)-xp(i))/(xp(i+1)-xp(i))*(pq(i+1)-pq(i));
        end
    end
    
%     % first order with matlab
%     v = interp1(xp,vq,x,'pchip');
%     p = interp1(xp,pq,x,'pchip');
         
    % other field variables
    for i=1:Nx
       epsz2(i) = 1/E2*(p(i)-p0)*Dn/2/dn*(2*epsz(i)+1)*(1-exp(-E2/eta2*dt_real)) + ...
           epsz2(i)*exp(-E2/eta2*dt_real);
       epsz(i) = epsz2(i) + 1/E1*(p(i)-p0)/2*Dn/dn*(2*epsz(i)+1)/(epsz(i)+1)^beta;
       a(i) = sqrt(E1*dn/rho/Dn*(epsz(i)+1)^beta);
       D(i) = Dn*(1+epsz(i));
       
    end
    
    if(mod(j,1) == 0)
        fprintf('Time: %6.4f s \n', t(j));
%             PLOT
        subplot(3,1,1) 
        grid on;
         plot(x,v,'linewidth',1.5);
        xlabel('x, mm');
        ylabel('Velocity, m/s');
        title(['time = ', num2str(t(j)), ' s'])

        subplot(3,1,2)
        grid on;
        plot(x,(p-p0)/mmHg_Pa,'linewidth',1.5);
        xlabel('x, mm');
        ylabel('Pressure, mmHg');

        subplot(3,1,3)
        grid on;
        plot(x,D/1000.,'linewidth',1.5);
        xlabel('x, mm');
        ylabel('Diameter, mm');
%         pause; 
    end 
    
    j=j+1;
    pt(:,j) = p;
    vt(:,j) = v;
    at(:,j) = a;
    Dt(:,j) = D;
    epszt(:,j) = epsz;
    epsz2t(:,j) = epsz2;
end

pt = pt(:,1:j);
vt = vt(:,1:j);
at = at(:,1:j);
Dt = Dt(:,1:j);
epszt = epszt(:,1:j);
epsz2t = epsz2t(:,1:j);









