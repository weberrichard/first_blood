clear;

par.l = 0.2;
par.ds = 0.011;
par.de = 0.01;
par.As = par.ds*par.ds*pi/4;
par.Ae = par.de*par.de*pi/4;

nx = 10;
dx = par.l/(nx-1);
x = linspace(0,par.l,nx);

An = par.As + x/par.l*(par.Ae-par.As);
dAn = 1/par.l*(par.Ae-par.As);

A = zeros(nx,1);
A(nx) = par.Ae;

for i=nx-1:-1:1
    A(i) = A(i+1) - dx*( sqrt(A(i+1)/An(i+1))*dAn + 2*sqrt(A(i+1))/An(i+1)*(sqrt(A(i+1))-sqrt(An(i+1))));
end

figure(2);
plot(x,A,'x');
hold on;
plot(x,An,'k');
