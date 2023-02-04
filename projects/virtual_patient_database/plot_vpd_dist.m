clear;

data = importdata('ki10000.csv');
data = data(2:end,2:end);

lit=[
65    , 12.84 % radial dias [mmHg]
123   , 22.75 % radial sys [mmHg]
75.6  , 12.7  % aortic dias [mmHg]
113   , 11.2  % aortic sys [mmHg]
75.58 , 6.01  % carotid dias [mmHg]
122.78, 12.01 % carotid sys [mmHg]
4570  , 1090  % cardiac output [ml/min]
7.63  , 1.563 % aortic PWV [m/s]
8.1   , 1.8   % car-fem PWV [m/s]
10.43 , 1.661 % bra-rad PWV [m/s]
9.794 , 1.778 % fem-ank PWV [m/s]
];

tit = ["Radial diastole", "Radial systole","Aortic diastole", "Aortic systole","Carotid diastole", "Carotid systole","Cardiac output","Aortic PWV", "Carotid-femoral PWV","Brachial radial PWV", "Femoral-ankle PWV"]; 

n = 10;

fig1 = figure('position',[100,200,1000,1000]);

% plotting pressures
for i=1:6
    % simulation data
    subplot(7,2,i);
    h1=histogram(data(i,:),n);
    h1.Normalization = 'pdf';
    hold on; grid on;
    title(tit(i));
    xlabel('Pressure [mmHg]');
    
    % literature data
    x = linspace(lit(i,1)-3*lit(i,2),lit(i,1)+3*lit(i,2),1000);
    y = pdf('Normal',x,lit(i,1),lit(i,2));
    plot(x,y,'linewidth',2);
end

% plotting cardiac output
% simulation data
subplot(7,2,[7 10]);
h1=histogram(data(8,:),n);
h1.Normalization = 'pdf';
hold on; grid on;
title(tit(7));
xlabel('Flow rate [ml/min]');

% literature data
x = linspace(lit(7,1)-3*lit(7,2),lit(7,1)+3*lit(7,2),1000);
y = pdf('Normal',x,lit(7,1),lit(7,2));
plot(x,y,'linewidth',2);


% pulse wave velocities
for i=1:4
    % simulation data
    subplot(7,2,i+10);
    h1=histogram(data(i+12,:),n);
    h1.Normalization = 'pdf';
    hold on; grid on;
    title(tit(i+7));
    xlabel('PWV [m/s]');
    
    % literature data
    x = linspace(lit(i+7,1)-3*lit(i+7,2),lit(i+7,1)+3*lit(i+7,2),1000);
    y = pdf('Normal',x,lit(i+7,1),lit(i+7,2));
    plot(x,y,'linewidth',2);
end

saveas(fig1,'plots/vpd_dist.fig','fig');
saveas(fig1,'plots/vpd_dist.eps','epsc');
saveas(fig1,'plots/vpd_dist.png','png');
