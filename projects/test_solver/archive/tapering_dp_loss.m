clear;
mu=0.0034;
v0=1;
rs=0.01;
beta=1*pi/180;
l=0.2;
re=rs+2*beta*l;
rho=1055;

n=100;

p0 = 1e5+100*133;

r=linspace(rs,re,n);
dx=l/(n-1);
x=linspace(0,l,n);

dr=(re-rs)/l;

v=zeros(1,n);
v(1) = v0;
dp1=zeros(1,n);
dp2=zeros(1,n);
dp3=zeros(1,n);
p1=zeros(1,n);
p2=zeros(1,n);
p3=zeros(1,n);
p1(1) = p0;
p2(1) = p0;
p3(1) = p0;
for i=1:n-1
    v(i+1) = rs^2*v0/(rs+x(i)/l*(re-rs))^2;
    dp1(i+1) = dp1(i) - dx*(4*mu*v(i+1)/r(i)^2 + 2*p1(i)/r(i)*dr + 2*rho*v(i+1)^2/r(i)*dr - 2*rho*v(i+1)*(v(i+1)-v(i))/dx);
    dp2(i+1) = dp2(i) - dx*(4*mu*v(i+1)/r(i)^2 + 2*p2(i)/r(i)*dr);
    dp3(i+1) = dp3(i) - dx*4*mu*v(i+1)/r(i)^2;
    
    p1(i+1)  = p0+rho/2*v0^2-rho/2*v(i+1)^2-dp1(i+1);
    p2(i+1)  = p0+rho/2*v0^2-rho/2*v(i+1)^2-dp2(i+1);
    p3(i+1)  = p0+rho/2*v0^2-rho/2*v(i+1)^2-dp3(i+1);
end

plot(x,dp1,'o')
hold on;
grid on;
plot(x,dp2,'x')
plot(x,dp3,'x')
legend(["all","only pres","original"]);

figure;
plot(x,p1,'o')
hold on;
grid on;
plot(x,p2,'x')
plot(x,p3,'x')
legend(["all","only pres","original"]);