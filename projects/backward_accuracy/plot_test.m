clear;

edge_id = 'A04';

A2 = readmatrix(['results/P045/',edge_id,'.txt']);
to = A2(:,1);
po = A2(:,2);
vo = A2(:,4);
do = A2(:,14);

% A2 = readmatrix('back_p_in.txt');
% tin = A2(:,1);
% pin = A2(:,2);
% 
% A2 = readmatrix('back_p_in_ip.txt');
% tip = A2(:,1);
% pip = A2(:,2);

B = readmatrix(['results/P045_back/',edge_id,'.txt']);
t = B(:,1);
p = B(:,2);
v = B(:,4);

po = interp1(to,po,t);

fig1=figure();
hold on;
% plot(tb,(pb-1e5)/133.3616);
plot(t,(p-1e5)/133.3616);
plot(t,(po-1e5)/133.3616);
% plot(t,p);
grid on;
xlim([0,5]);
xlabel('time [s]');
ylabel('pressure [mmHg]');
legend('forward','backward');
% saveas(fig1,'pressure.png','png');

e = (p-po)./po;
fig2=figure();
plot(t,e*100);
grid on;
xlim([0,5]);
xlabel('time [s]');
ylabel('error [%]');
% saveas(fig2,'error.png','png');

% figure();
% plot(t,v);
% hold on;
% plot(ti,vi);
% grid on;

% % figure();
% plot(t,d);
% hold on;
% plot(ti,di);
% grid on;