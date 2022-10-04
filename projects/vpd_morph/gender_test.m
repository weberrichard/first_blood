clear;

N = 1000;
pm = 110;
ps = 30;
p = ps*randn(N,1)+pm;

pm1 = 115;
ps1 = 30;

pm2 = 105;
ps2 = 30;

gend = randi(2,N,1);

r = round(log(N)/log(2)+1);
b = prctile(p,linspace(0,100,r));

% orig xhi
p1 = p(gend==1);
p2 = p(gend==2);
N1 = length(p1);
N2 = length(p2);

A1 = p1>b;
nu1 = A1(:,1:end-1)-A1(:,2:end);
nu1(1) = nu1(1)+1;
nu1 = sum(nu1);
pi1 = normcdf(b,pm1,ps1);
pi1(1) = 0;
pi1(end) = 1;
pi1 = pi1(2:end) - pi1(1:end-1);
xhi1 = sum((nu1-N1*pi1).^2./N1./pi1);

A2 = p2>b;
nu2 = A2(:,1:end-1)-A2(:,2:end);
nu2(1) = nu2(1)+1;
nu2 = sum(nu2);
pi2 = normcdf(b,pm2,ps2);
pi2(1) = 0;
pi2(end) = 1;
pi2 = pi2(2:end) - pi2(1:end-1);
xhi2 = sum((nu2-N2*pi2).^2./N2./pi2);

xhi_orig = xhi1 + xhi2;

good_change = zeros(N,1);

for i=1:N
    % changing gender
    gend_orig = gend(i);
    if(gend(i)==1)
        gend(i) = 2;
    else
        gend(i) = 1;
    end

    p1 = p(gend==1);
    p2 = p(gend==2);
    N1 = length(p1);
    N2 = length(p2);

    A1 = p1>b;
    nu1 = A1(:,1:end-1)-A1(:,2:end);
    nu1(1) = nu1(1)+1;
    nu1 = sum(nu1);
    pi1 = normcdf(b,pm1,ps1);
    pi1(1) = 0;
    pi1(end) = 1;
    pi1 = pi1(2:end) - pi1(1:end-1);
    xhi1 = sum((nu1-N1*pi1).^2./N1./pi1);
    
    A2 = p2>b;
    nu2 = A2(:,1:end-1)-A2(:,2:end);
    nu2(1) = nu2(1)+1;
    nu2 = sum(nu2);
    pi2 = normcdf(b,pm2,ps2);
    pi2(1) = 0;
    pi2(end) = 1;
    pi2 = pi2(2:end) - pi2(1:end-1);
    xhi2 = sum((nu2-N2*pi2).^2./N2./pi2);
    
    xhi = xhi1 + xhi2;

    if(xhi<xhi_orig)
        good_change(i) = 1;
    end

    % changing back
    gend(i) = gend_orig;
end

for i=1:N
    if(good_change(i)==1)
        if(gend(i)==1)
            gend(i) = 2;
        else
            gend(i) = 1;
        end
    end
end
p1 = p(gend==1);
p2 = p(gend==2);
N1 = length(p1);
N2 = length(p2);

A1 = p1>b;
nu1 = A1(:,1:end-1)-A1(:,2:end);
nu1(1) = nu1(1)+1;
nu1 = sum(nu1);
pi1 = normcdf(b,pm1,ps1);
pi1(1) = 0;
pi1(end) = 1;
pi1 = pi1(2:end) - pi1(1:end-1);
xhi1 = sum((nu1-N1*pi1).^2./N1./pi1);

A2 = p2>b;
nu2 = A2(:,1:end-1)-A2(:,2:end);
nu2(1) = nu2(1)+1;
nu2 = sum(nu2);
pi2 = normcdf(b,pm2,ps2);
pi2(1) = 0;
pi2(end) = 1;
pi2 = pi2(2:end) - pi2(1:end-1);
xhi2 = sum((nu2-N2*pi2).^2./N2./pi2);

xhi_final = xhi1 + xhi2;

g = repmat({'First'},N,1);
g1 = repmat({'Second'},N1,1);
g2 = repmat({'Third'},N2,1);
g = [g; g1; g2];

pp = [p;p1;p2];

figure;
boxplot(pp,g);
hold on;


