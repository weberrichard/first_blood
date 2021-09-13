clear;
% main parameters
HR  = 60;  % heart rate
per = 5; % number of simulation periods
mmHg_Pa = 133.322; % for converts

% creating the elastance function
T0 = 2.5; 
nt = 1000; % number of timesteps in a period
tn = linspace(0,T0,nt);
En = 1.55*((tn/0.7).^1.9./(1+(tn/0.7).^1.9)) .* (1./(1+(tn/1.17).^21.9));

% figure;
% plot(tn,En);
% grid on;

% real elastance
Emin = 0.06; % max elastance
Emax = 2.5; % min elastance

E0 = (Emax-Emin)*En + Emin;
t0 = linspace(0,60/HR,nt);
%making periodic
E=[];
t=[];
te=t0(end);
for i=1:per
   E = [E,E0];
   if(i==1)
       t=t0;
   else
       t = [t,t0+max(t)];
   end
end
E=[E,E(1)];
t=[t,max(t)+t(2)];

% convert E
E = E * mmHg_Pa * 1000 * 1000;

% E=E/100;

% figure;
% plot(t,E);
% grid on;

% parameters
R1 = 1.49955; % mmHg s / ml   0.5-2.0
R2 = 0.00499852; % mmHg s / ml
R3 = 0.001; % 0.001 mmHg s / ml
R4 = 0.0398; % mmHg s / ml
% C1 = 1/E(1); % ml / mmHg
C2 = 4.40133; % ml / mmHg
C3 = 1.33039; % ml / mmHg
L = 0.0005; % mmHg s2 / ml
V0 = 10; % cm3

% converting to SI
R1 = R1 * mmHg_Pa*1000*1000; % Pa/m3
R2 = R2 * mmHg_Pa*1000*1000; % Pa/m3
R3 = R3 * mmHg_Pa*1000*1000; % Pa/m3
R4 = R4 * mmHg_Pa*1000*1000; % Pa/m3
% C1 = C1; % m3/Pa
C2 = C2 / mmHg_Pa/1000/1000; % m3/Pa
C3 = C3 / mmHg_Pa/1000/1000; % m3/Pa
L = L * mmHg_Pa*1000*1000; % Pas2/m3
V0 = V0 /1000/1000; % m3

% initial conditions
LVV0 = 61; % ml
Pla0 = 6.318; %6.318; %6.318; % mmHg
% Pa0 = 46; %45.891 ; % mmHg
qa0 = 0; % m3/s
Plv0 = E(1)*(LVV0/1000/1000-V0)/mmHg_Pa; % mmHg
Pa0 = 90; %45.891 ; % mmHg

LVV = zeros(1,nt*per+1); LVV(1) = LVV0;
Pla = zeros(1,nt*per+1); Pla(1) = Pla0;
Pa = zeros(1,nt*per+1); Pa(1) = Pa0;
qa = zeros(1,nt*per+1); qa(1) = qa0;
Plv = zeros(1,nt*per+1); Plv(1) = Plv0;
EC1 = zeros(1,nt*per+1);
EC2 = zeros(1,nt*per+1);
EC3 = zeros(1,nt*per+1);
EL  = zeros(1,nt*per+1);
Ener = zeros(1,nt*per+1);

z = zeros(4,nt*per+1);
z(:,1) = [LVV0/1000/1000-V0; Pla0* mmHg_Pa; Pa0* mmHg_Pa; qa0];

% energy of the system in SI
EC1(1) = .5*E(1)*z(1,1)^2;
EC2(1) = .5*C2*z(2,1)^2;
EC3(1) = .5*C3*z(3,1)^2;
EL(1)  = .5*L*z(4,1)^2;
Ener(1) = EC1(1) + EC2(1) + EC3(1) + EL(1);

% ejection
Ae = [0,0,0,-1;0,-1/R1/C2,1/R1/C2,0;0,1/R1/C3,-1/R1/C3,1/C3;E(1)/L,0,-1/L,-(R3+R4)/L];

% filling
Af = [-E(1)/R2,1/R2,0,0;E(1)/R2/C2,-(R1+R2)/R1/R2/C2,1/R1/C2,0;0,1/R1/C3,-1/R1/C3,0;0,0,0,0];

% relaxation
Ar = [0,0,0,0;0,-1/R1/C2,1/R1/C2,0;0,1/R1/C3,-1/R1/C3,0;0,0,0,0];

fig1=figure;
hold on; grid on;
% calculation
for i=1:nt*per
    
    % ejection
%     if(Plv(i)>Pa(i) && Plv(i)>Pla(i))
%         z(:,i+1) = z(:,i) + (t(i+1)-t(i))*Ae*z(:,i);
% %         disp('eject');
%     end
    
    % filling
%     if(Plv(i)<Pa(i) && Plv(i)<Pla(i))
        z(:,i+1) = z(:,i) + (t(i+1)-t(i))*Af*z(:,i);
%          disp('fill');
%     end
     
    % relaxation
%     if(Plv(i)<Pa(i) && Plv(i)>Pla(i))
%         z(:,i+1) = z(:,i) + (t(i+1)-t(i))*Ar*z(:,i);
% %         disp('relax');
%     end
%     
    % matrix update
    Ae(4,1) = E(i+1)/L;
    Af(1,1) = -E(i+1)/R2;
    Af(2,1) = E(i+1)/R2/C2;
    
    % field var save
    LVV(i+1) = (z(1,i+1)+V0)*1000*1000; % ml
    Pla(i+1) = z(2,i+1)/mmHg_Pa; % mmHg
    Pa(i+1) = z(3,i+1)/mmHg_Pa; % mmHg
    qa(i+1) = z(4,i+1)*1000*1000; % ml
    Plv(i+1) = E(i+1)*z(1,i+1)/mmHg_Pa; % mmHg
    
    % energy of the system in SI
    EC1(i+1) = .5*E(i+1)*z(1,i+1)^2;
    EC2(i+1) = .5*C2*z(2,i+1)^2;
    EC3(i+1) = .5*C3*z(3,i+1)^2;
    EL(i+1)  = .5*L*z(4,i+1)^2;
    Ener(i+1) = EC1(i+1) + EC2(i+1) + EC3(i+1) + EL(i+1);
    
    % PLOT
%     if(mod(i,100 )==0)
%         cla;
%         plot(t(1:i),Pla(1:i));
%         grid on; hold on; 
%         plot(t(1:i),Plv(1:i));
%         plot(t(1:i),Pa(1:i));
%         legend('Pla','Plv','Pa');
%         pause;
%     end
end

% cla;
plot(t(1:i),Pla(1:i),'linewidth',1.5);
grid on; hold on; 
plot(t(1:i),Plv(1:i),'linewidth',1.5);
plot(t(1:i),Pa(1:i),'linewidth',1.5);
% plot(t(1:i),LVV(1:i),'linewidth',1.5);
legend('Pla','Plv','Pa','location','northeast');
xlabel('time [s]');
ylabel('Pressure [mmHg], LVV [ml]');

% fig2=figure;
% plot(LVV(end-nt:end),Plv(end-nt:end),'linewidth',1.5);
% grid on
% xlabel('LVV [ml]');
% ylabel('Plv [mmHg]');

% saveas(fig1,'time-pressure.png','png');
% saveas(fig2,'lvv-plv.png','png');

% energy of elastance capacitor
% E = E / mmHg_Pa / 1e6;
% Ec = .5*1./E.*Plv.^2;

% fig3=figure;
% plot(t,Ec);
% % plot(t,1./E);
% % yyaxis right
% % plot(t,Plv);
% grid on;


% fig4 = figure;
% plot(t,Ener,'linewidth',1.5);
% grid on
% hold on
% plot(t,EC1,'linewidth',1.5);
% plot(t,EC2,'linewidth',1.5);
% plot(t,EC3,'linewidth',1.5);
% plot(t,EL,'linewidth',1.5);
% legend('total','C1','C2','C3','L');



