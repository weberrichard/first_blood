clear; close all;

file_name = 'log_de_5.txt';
d = importdata(file_name);

literature = [65,123,65,103,75.58,122.78,350.4,4570,235,385,75,142.8,7.63,8.1,10.43,9.79];
dev = [12.84,22.75,7.14,8.24,6.01,12.01,126.6,1090,94.328,99,40.2,70.8,1.56,1.80,1.66,1.78];

best_n = length(d);

ff = d(:,1);
[~,idx] = sort(ff);
ff = ff(idx(1:best_n));

% outputs
p = d(idx(1:best_n),2:7);
q = d(idx(1:best_n),8:13);
a = d(idx(1:best_n),14:17);
out = d(idx(1:best_n),2:17);
% inputs, parameters
heart = d(idx(1:best_n),18:20);
res = d(idx(1:best_n),21:32);
cap = d(idx(1:best_n),33:36);
len = d(idx(1:best_n),37:40);
dia = d(idx(1:best_n),41:44);
thi = d(idx(1:best_n),45:48);
ela = d(idx(1:best_n),49:56);
par = d(idx(1:best_n),18:56);

tc_perif = res(:,5:8).*cap;

figure('position',[0,50,300,600]);
boxplot(ff);
title('Fitness function [-]');

figure('position',[350,50,1200,600]);
subplot(1,5,[1,2]);
boxplot(p);
hold on;
plot([linspace(.5,5.5,6)',linspace(1.5,6.5,6)']',[p(1,:)',p(1,:)']','m','linewidth',1.5);
plot([linspace(.5,5.5,6)',linspace(1.5,6.5,6)']',[literature(1:6)',literature(1:6)']','g','linewidth',1.5);
plot([linspace(.5,5.5,6)',linspace(1.5,6.5,6)']',[literature(1:6)'+dev(1:6)',literature(1:6)'+dev(1:6)']','--g','linewidth',1.2);
plot([linspace(.5,5.5,6)',linspace(1.5,6.5,6)']',[literature(1:6)'-dev(1:6)',literature(1:6)'-dev(1:6)']','--g','linewidth',1.2);
ylim([50,150]);
title('Pressure [mmHg]');

subplot(1,5,3);
boxplot(q(:,[1,3:end]));
hold on;
plot([linspace(.5,4.5,5)',linspace(1.5,5.5,5)']',[q(1,[1,3:end])',q(1,[1,3:end])']','m','linewidth',1.5);
plot([linspace(.5,4.5,5)',linspace(1.5,5.5,5)']',[literature([7,9:12])',literature([7,9:12])']','g','linewidth',1.5);
plot([linspace(.5,4.5,5)',linspace(1.5,5.5,5)']',[literature([7,9:12])'+dev([7,9:12])',literature([7,9:12])'+dev([7,9:12])']','--g','linewidth',1.2);
plot([linspace(.5,4.5,5)',linspace(1.5,5.5,5)']',[literature([7,9:12])'-dev([7,9:12])',literature([7,9:12])'-dev([7,9:12])']','--g','linewidth',1.2);
title('Volume flow rate [ml/min]');

subplot(1,5,4);
boxplot(q(:,2));
hold on;
plot([linspace(.5,0.5,1)',linspace(1.5,1.5,1)']',[q(1,2)',q(1,2)']','m','linewidth',1.5);
plot([linspace(.5,0.5,1)',linspace(1.5,1.5,1)']',[literature(8)',literature(8)']','g','linewidth',1.5);
plot([linspace(.5,0.5,1)',linspace(1.5,1.5,1)']',[literature(8)'+dev(8)',literature(8)'+dev(8)']','--g','linewidth',1.2);
plot([linspace(.5,0.5,1)',linspace(1.5,1.5,1)']',[literature(8)'-dev(8)',literature(8)'-dev(8)']','--g','linewidth',1.2);
ylim([3000,6000]);
title('Cardiac output [ml/s]');

subplot(1,5,5);

boxplot(a);
hold on;
plot([linspace(.5,3.5,4)',linspace(1.5,4.5,4)']',[a(1,:)',a(1,:)']','m','linewidth',1.5);
plot([linspace(.5,3.5,4)',linspace(1.5,4.5,4)']',[literature(13:16)',literature(13:16)']','g','linewidth',1.5);
plot([linspace(.5,3.5,4)',linspace(1.5,4.5,4)']',[literature(13:16)'+dev(13:16)',literature(13:16)'+dev(13:16)']','--g','linewidth',1.2);
plot([linspace(.5,3.5,4)',linspace(1.5,4.5,4)']',[literature(13:16)'-dev(13:16)',literature(13:16)'-dev(13:16)']','--g','linewidth',1.2);
ylim([6,13]);
title('Pulse wave velocity [m/s]');

figure('position',[0,750,1500,600]);
subplot(2,5,1);
boxplot(heart);
hold on;
title('Heart: Rmit, Raor, Patr');
plot([linspace(.5,2.5,3)',linspace(1.5,3.5,3)']',[heart(1,:)',heart(1,:)']','m','linewidth',1.5);

subplot(2,5,[2:4]);
boxplot(res);
hold on;
title('Rnode, Rperif1, Rperif2');
plot([linspace(.5,11.5,12)',linspace(1.5,12.5,12)']',[res(1,:)',res(1,:)']','m','linewidth',1.5);

subplot(2,5,5);
boxplot(cap);
hold on;
title('Cperif');
plot([linspace(.5,3.5,4)',linspace(1.5,4.5,4)']',[cap(1,:)',cap(1,:)']','m','linewidth',1.5);

subplot(2,5,6);
boxplot(len);
hold on;
title('Length');
plot([linspace(.5,3.5,4)',linspace(1.5,4.5,4)']',[len(1,:)',len(1,:)']','m','linewidth',1.5);

subplot(2,5,7);
boxplot(dia);
hold on;
title('Diameter');
plot([linspace(.5,3.5,4)',linspace(1.5,4.5,4)']',[dia(1,:)',dia(1,:)']','m','linewidth',1.5);

subplot(2,5,8);
boxplot(thi);
hold on;
title('Thickness');
plot([linspace(.5,3.5,4)',linspace(1.5,4.5,4)']',[thi(1,:)',thi(1,:)']','m','linewidth',1.5);

subplot(2,5,[9,10]);
boxplot(ela);
hold on;
title('Espring, Evoigt');
plot([linspace(.5,7.5,8)',linspace(1.5,8.5,8)']',[ela(1,:)',ela(1,:)']','m','linewidth',1.5);


