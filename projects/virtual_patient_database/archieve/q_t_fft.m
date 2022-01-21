clear; close all;

v = importdata("q-t_diagram.txt");
to = v(:,1)/0.8;
xo = v(:,2);

figure;
plot(to,xo,'x');
hold on
grid on

nt = 1000;
t = 0:1/nt:1;
x = interp1(to,xo,t,'PCHIP');
plot(t,x);
legend('original','resampled');

x = x(1:401);
t = t(1:401);

N = 5;
xN = repmat(x(1:end-1),1,N);
xN = [xN,x(end)];
tN = linspace(0,2,2001);

f2 = fit(tN',xN','fourier3')
plot(f2);
coef = coeffvalues(f2);
a0 = coef(1);
a = coef(2:2:end-1);
b = coef(3:2:end);
cint = confint(f2);
a0b = cint(:,1);
ab = cint(:,2:2:end-1);
bb = cint(:,3:2:end);
w = f2.w;
xf = zeros(size(t));
xh = xf;
xl = xf;

xf = xf + a0;
xl = xl + a0b(1);
xh = xh + a0b(2);
for i=1:length(a)
    xf = xf + a(i)*cos(i*t*w) + b(i)*sin(i*t*w);
    xl = xl + ab(1,i)*cos(i*t*w) + bb(1,i)*sin(i*t*w);
    xh = xh + ab(2,i)*cos(i*t*w) + bb(2,i)*sin(i*t*w);
end
plot(t,xf,'linewidth',2);
plot(t,xl,'k');
plot(t,xh,'k');

% creating q-t function
t = 0:.001:1;
q = zeros(size(t));

delta_t = 0.0059;
t = t-delta_t;
q = q + coef(1);
for i=1:length(a)
    q = q + coef(2*i)*cos(i*t*w) + coef(2*i+1)*sin(i*t*w);
end
t = t+delta_t;

q(t>.4) = 0;

figure;
grid on; hold on;
plot(t,q,'linewidth',2);
xlabel("time [-]");
ylabel("volume flow rate [ml/s]");

% q with noise
n = 5;
sigma = .2;

for k=1:n
   q = zeros(size(t));
   noise = randn(size(coef))*sigma+1;
   coefn = coef.*noise;
   q = q + coefn(1);
   t = t-delta_t;
   for i=1:length(a)
    q = q + coef(2*i)*cos(i*t*w) + coef(2*i+1)*sin(i*t*w);
   end
   t = t+delta_t;
   q(t>.4) = 0;
   plot(t,q);
end
    

    
    
    
    
