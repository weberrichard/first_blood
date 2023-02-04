% PLAN B

rho = 1055;
Q0 = 15*1e-6;
A1 = 1.98857e-6;
A2 = 1.85632e-5;

d1 = sqrt(4*A1/pi);
d2 = sqrt(4*A2/pi);

Q = (0:.1:20)*1e-6;

dpv = 8*rho*Q.^2/pi^2*(1/d1^2 - 1/d2^2)^2;

R = 16*rho/pi^2*(1./d1.^2 - 1./d2.^2)^2*Q0;
R2 = 8*rho/pi^2*(1./d1.^2 - 1./d2.^2)^2;

dp0 = 8*rho*Q0^2/pi^2*(1/d1^2 - 1/d2^2)^2;
Q_lin = (5:.1:15)*1e-6;
dp_lin = dp0 + R*(Q_lin-Q0);

figure;
plot(Q,dpv);
hold on; grid on;
plot(Q_lin,dp_lin);
plot(Q,R2*Q.^2,'x');