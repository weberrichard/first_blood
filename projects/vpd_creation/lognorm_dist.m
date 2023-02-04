clear;

sigma = 0.472;
mu = -sigma^2/2; % to get 1 as expected value

n = 1e5;

v = lognrnd(mu,sigma,n,1);

M = exp(mu+sigma^2/2);
S = sqrt((exp(sigma^2)-1)*exp(2*mu+sigma^2));


sv = 0.1:0.001:2;
mv = -sv.^2./2;
ss = sqrt((exp(sv.^2)-1).*exp(2*mv+sv.^2));

figure;
plot(sv,ss,'linewidth',2);
grid on;
xlabel('sigma');
ylabel('deviation');
