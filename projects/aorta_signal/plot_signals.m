clear; close all;

files = dir('calculated_aorta_signals');

% plotting original
figure;
hold on; grid on;
for i=3:length(files)
    data = importdata(['calculated_aorta_signals/',files(i).name]);
    t = data(:,1);
    p = data(:,2);
    plot(t,p);
end

saveas(gca,'plots/orig.png','png');

n = 1001;
tt = linspace(0,1,n);
pp = zeros(length(files)-2,n);

% normalizing stuff
figure;
hold on; grid on;
for i=3:length(files)
    data = importdata(['calculated_aorta_signals/',files(i).name]);
    t = data(:,1);
    p = data(:,2);
    t = t./max(t);
    p = (p-min(p))./(max(p)-min(p));
    pp(i-2,:) = interp1(t,p,tt);
    plot(t,p,'r');
end

plot(tt,mean(pp),'k','linewidth',3);

saveas(gca,'plots/normalized.png','png');

mpp = mean(pp);
v = [tt',mpp'];

