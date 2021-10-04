clear;

tn = 0:.001:1;
En = 17.4073 * tn.^1.9 ./ (1+11.2305*tn.^1.9) * 1 ./ (1+1.6658e7*tn.^21.9);

figure;
plot(tn,En,'linewidth',2);
grid on;

Emin = 0.06;
Emax = 2.5;
HR = 75;

t = linspace(0,60/HR,length(tn));
E = (Emax-Emin)*En + Emin;

figure;
plot(t,E,'linewidth',2);
grid on;