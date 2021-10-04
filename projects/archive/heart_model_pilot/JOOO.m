clear;

n=4;
m=5;
e=1;

A = zeros(n+m+2*e);
b = zeros(n+m+2*e,1);

Tmax = 5;
dt = 1e-3;
nt = Tmax/dt;
t = zeros(nt,1);
i=1;

R1 = 1.5;
R2 = 5e-3;
C2 = 4.4;
C3 = 1.33;

% N1, N2, N3, G
p = [6.318;3.06;90;0];
% R1, R2, C2, C3, 1/E
q = [0;0;0;0;0];
% y1,y2
E0 = elastance(0);
y = [p(4)/E0;p(2)/E0];

pN1 = zeros(nt,1); pN1(1) = p(1);
pN2 = zeros(nt,1); pN2(1) = p(2);
pN3 = zeros(nt,1); pN3(1) = p(3);
E = zeros(nt,1);

%    1    2   3   4   5  6   7   8  9  10 11
% x=[qR1,qR2,qC2,qC3,qE,pN1,pN2,pN3,pG,y1,y2]

while(t<Tmax)
    A(1,1) = R1;
    A(1,6) = -1; %N1
    A(1,8) = 1; %N3
    
    A(2,2) = R2;
    A(2,6) = -1; %N1
    A(2,7) = 1; %N2
    
    A(3,3) = dt/C2;
    A(3,6) = 1; %N1
    A(3,9) = -1; %G 
    
    A(4,4) = dt/C3;
    A(4,8) = 1; %N3
    A(4,9) = -1; %G
    
    A(5,5) = dt;
    A(5,10) = -1; % y1
    A(5,11) = 1; % y2
    
    A(6,1) = -1; %R1
    A(6,2) = -1; %R2
    A(6,3) = 1; %C2
    
    A(7,2) = 1; %R2
    A(7,5) = 1; %1/E
    
    A(8,1) = 1; %R1
    A(8,4) = 1; % C3
    
    A(9,9) = 1; %G
        
    E(i) = elastance(t(i));
    A(10,10) = E(i);
    A(10,9) = -1;
    
    A(11,11) = E(i);
    A(11,7) = -1;
    
    b(3) = p(1)-p(4);
    b(4) = p(3)-p(4);
    b(5) = y(2)-y(1);
    
    x = A\b;
    
    q = x(1:5);
    p = x(6:9);
    y = x(10:11);
    
    pN1(i+1) = p(1);
    pN2(i+1) = p(2);
    pN3(i+1) = p(3);
    
    t(i+1) = t(i)+dt;
    i = i+1;
end

figure;
plot(t,pN1,'linewidth',1.5);
hold on; grid on;
plot(t,pN2,'linewidth',1.5);
plot(t,pN3,'linewidth',1.5);
legend('N1','N2','N3');

function out=elastance(t)

    heart_rate = 60;
    elastance_max = 2.5;
    elastance_min = 0.06;

	tn = t * heart_rate/60;

	while(tn>1)
		tn = tn - 1.;
	end

	En = 17.4073 * tn^1.9 / (1.+11.2305*tn^1.9) * 1. / (1.+1.6658e7*tn^21.9);

	E = (elastance_max-elastance_min)*En + elastance_min;
    
    out = E;
end


