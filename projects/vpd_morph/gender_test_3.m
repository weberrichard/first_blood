clear;

N2 = 100;
N1 = 100;
N = N1+N2;
n = 200;
Nit = 20;

% pm = 100;
% ps = 50;
% p = ps*randn(N,1)+pm;

pm1 = 120;
ps1 = 30;

pm2 = 80;
ps2 = 30;

p1o = ps1*randn(N1,1)+pm1;
p2o = ps2*randn(N2,1)+pm2;
p=[p1o;p2o];

% adressing some genders
gend = zeros(N,1);
gend(1:n) = randi(2,n,1);

% orig xhi
p1 = p(gend==1);
p2 = p(gend==2);

xhi1 = xhi_test(p1,pm1,ps1);
xhi2 = xhi_test(p2,pm2,ps2);
% xhi1 = mean_std_dif(p1,pm1,ps1);
% xhi2 = mean_std_dif(p2,pm2,ps2);
xhi_orig = xhi1 + xhi2;

good_change = zeros(N,1);

% adressing more sexes
for j=1:Nit
    for i=1:N
        if(j>1 || i>n)
            % changing gender
            gend_orig = gend(i);
            if(gend(i)==1)
                gend(i) = 2;
            else
                gend(i) = 1;
            end
            
            p1 = p(gend==1);
            p2 = p(gend==2);
            xhi11 = xhi_test(p1,pm1,ps1);
            xhi12 = xhi_test(p2,pm2,ps2);
    %         xhi11 = mean_std_dif(p1,pm1,ps1); 
    %         xhi12 = mean_std_dif(p2,pm2,ps2);
            xhi = xhi11 + xhi12;
            
            if(xhi<xhi_orig)
                good_change(i) = 1;
            end

            % changing back
            gend(i) = gend_orig;
        end
    end
    % making the change actually
    for i=1:N
        if(good_change(i)==1)
            if(gend(i)==1)
                gend(i) = 2;
            else
                gend(i) = 1;
            end
        end
    end
end

% final evaluation
p1 = p(gend==1);
p2 = p(gend==2);
N1 = length(p1);
N2 = length(p2);

xhi_final = xhi_test(p1,pm1,ps1) + xhi_test(p2,pm2,ps2);
% xhi_final = mean_std_dif(p1,pm1,ps1) + mean_std_dif(p2,pm2,ps2);

g = repmat({'First'},N,1);
g1 = repmat({'Second'},N1,1);
g2 = repmat({'Third'},N2,1);
g = [g; g1; g2];

pp = [p;p1;p2];

figure;
boxplot(pp,g);
hold on; grid on;

function out=xhi_test(p,m,s)
    N = length(p);
    r = round(log(N)/log(2)+1);
    b = prctile(p,linspace(0,100,r));
    A = p>b;
    nu = A(:,1:end-1)-A(:,2:end);
    nu(1) = nu(1)+1;
    nu = sum(nu);
    pi = normcdf(b,m,s);
    pi(1) = 0;
    pi(end) = 1;
    pi = pi(2:end) - pi(1:end-1);
    out = sum((nu-N*pi).^2./N./pi);
end

function out=mean_std_dif(p,m,s)
    out = abs(mean(p)-m)/m + abs(std(p)-s)/s;
end
