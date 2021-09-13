clear;
T0 = 2.5; 
nt = 2500; % number of timesteps in a period
tn = linspace(0,T0,nt);
En = 1.55*((tn/0.7).^1.9./(1+(tn/0.7).^1.9)) .* (1./(1+(tn/1.17).^21.9));

Emin = 0.06; % max elastance
Emax = 2.5; % min elastance
En = (Emax-Emin)*En + Emin;

En = [En,En,En,En,En];
tn = [tn,tn+T0+T0/nt,tn+2*T0,tn+3*T0,tn+4*T0];
t = tn/2.5;

% derivative
Ep = diff(En)/(tn(2)-tn(1));
EpE = Ep./En(1:end-1);

Ep_anal = (1449.81*tn.^0.9 - 4.69371165466148e-13*tn.^2.8 - 490.138*tn.^22.8 - 1056.93*tn.^24.7)./((0.507792 + tn.^1.9).^2.*(31.1365 + tn.^21.9).^2);
Ep_anal = (Emax-Emin)*Ep_anal;

figure;
plot(t,En,'linewidth',3);
grid on; hold on;
% yyaxis right
plot(t(1:end-1),Ep,'linewidth',3);
% plot(tn(1:end-1),EpE,'linewidth',3);
% plot(tn(1:end-1),EpE,'linewidth',3);
plot(t(1:end-1),EpE,'linewidth',3);
dEpE = diff(EpE)/(tn(2)-tn(1));
plot(t(2:end-1),dEpE,'linewidth',3);
% plot(tn,Ep_anal,'x');
% legend('E','Ep','Ep/E');
legend('E','Ep','Ep/E','(Ep/E)p');
% legend('Ep/E','(Ep/E)p');
xlim([0,1]);
xlabel('time [s]','fontsize',14);
ylabel('elastance [mmHg/ml]','fontsize',14);

% set(gca,'XScale','log','YScale','log');
