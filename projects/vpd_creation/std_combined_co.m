clear;

a=[5439.9 5462.2 5261.5 5038.6 4737.6 4447.8];
s=[652.8 652.8 612.4 589.2 548.8 496.8];
n=[1200,1200,1200,1200,1200,1200];

sk = sdc_gen(a,n,s)

function s1 = sdc_gen(x,n,s)
    x1 = x(1);
    n1 = n(1);
    s1 = s(1);
    for i=2:length(n)
        s2 = sdc(n1,n(i),s1,s(i),x1,x(i));
        x2 = (n1*x1+n(i)*x(i))/(n1+n(i));
        n2 = n1+n(i);

        s1 = s2;
        x1 = x2;
        n1 = n2;
    end
end

function s = sdc(n1,n2,s1,s2,m1,m2)
    s = sqrt(((n1-1)*s1^2+(n2-1)*s2^2)/(n1+n2-1)+n1*n2*(m1-m2)^2/(n1+n2)/(n1+n2-1));
end