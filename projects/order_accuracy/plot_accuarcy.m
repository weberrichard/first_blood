clear;

A=[1.08876e+05, 1.13910e+05, 8.92452e-05, 1.71635e+00
1.08895e+05, 1.13870e+05, 8.98408e-05, 4.38393e+00
1.08911e+05, 1.13854e+05, 9.00229e-05, 1.11027e+01
1.08920e+05, 1.13848e+05, 9.01136e-05, 3.18025e+01
1.08924e+05, 1.13845e+05, 9.01589e-05, 9.94007e+01];

mmHg_to_Pa = 133.3616;

p_dia = (A(:,1)-1e5)/mmHg_to_Pa;
p_sys = (A(:,2)-1e5)/mmHg_to_Pa;
q_car = A(:,3)*1e3*60;
t_run = A(:,4);

n = 1101;
nx = [n,2*n,4*n,8*n,16*n]';
nx2 = linspace(n,16*n);

f_dia = fit(nx,p_dia,'power2');
p_dia = abs(p_dia-f_dia.c)/f_dia.c*100;
f_dia.a = abs(f_dia.a)/f_dia.c*100;
f_dia.c = 0;

f_sys = fit(nx,p_sys,'power2');
p_sys = abs(p_sys-f_sys.c)/f_sys.c*100;
f_sys.a = abs(f_sys.a)/f_sys.c*100;
f_sys.c = 0;

f_car = fit(nx,q_car,'power2');
q_car = abs(q_car-f_car.c)/f_car.c*100;
f_car.a = abs(f_car.a)/f_car.c*100;
f_car.c = 0;

x_time = linspace(min(t_run),max(t_run),1000);

fig=figure('position',[400,200,1000,800]);
subplot(2,2,1);
grid on;
hold on;
plot(nx2,f_dia(nx2),'r','linewidth',1.5);
plot(nx,p_dia,'xb','linewidth',2,'markersize',10);
xlabel('Division points [pc.]');
ylabel('Relative error [%]');
title('Aortic diastolic pressure');
legend('Power law fit','Simulation');
set(gca,'fontsize',14);
% set(gca,'xscale','log','yscale','log');

subplot(2,2,2);
grid on;
hold on;
plot(nx2,f_sys(nx2),'r','linewidth',1.5);
plot(nx,p_sys,'xb','linewidth',2,'markersize',10);
xlabel('Division points [pc.]');
ylabel('Relative error [%]');
title('Aortic systolic pressure');
legend('Power law fit','Simulation');
set(gca,'fontsize',14);
% set(gca,'xscale','log','yscale','log');

subplot(2,2,3);
grid on;
hold on;
plot(nx2,f_car(nx2),'r','linewidth',1.5);
plot(nx,q_car,'xb','linewidth',2,'markersize',10);
xlabel('Division points [pc.]');
ylabel('Relative error [%]');
title('Cardiac output');
legend('Power law fit','Simulation');
set(gca,'fontsize',14);
% set(gca,'xscale','log','yscale','log');

subplot(2,2,4);
f = fit(nx,t_run,'poly2');
x = linspace(n,16*n,100);
y = f.p1*x.^2 + f.p2*x + f.p3;
hold on;
plot(x,y,'r','linewidth',1.5);
plot(nx,t_run,'+b','linewidth',2,'markersize',10);
grid on;
xlabel('Division points [pc.]');
ylabel('Run time per cycle [s]');
title('Run time');
legend('Quadratic fit','Simulation','location','nw');
set(gca,'fontsize',14);

saveas(fig,'plots/order_accuracy.eps','epsc');
saveas(fig,'plots/order_accuracy.png','png');
saveas(fig,'plots/order_accuracy.fig','fig');