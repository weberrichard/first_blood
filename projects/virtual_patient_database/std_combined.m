clear;

n1 = 3837;
n2 = 1588;
n3 = 1588;

s1 = 1.7;
s2 = 1.363;
s3 = 1.303;

m1 = 10.7;
m2 = 9.78;
m3 = 8.57;

s12 = sdc(n1,n2,s1,s2,m1,m2);
m12 = (n1*m1+n2*m2)/(n1+n2);
n12 = n1+n2;
sc = sdc(n12,n3,s12,s3,m12,m3);
mc = (n12*m12+n3*m3)/(n12+n3);

ss = sdc_gen([m1,m2,m3],[n1,n2,n3],[s1,s2,s3]);

function s = sdc(n1,n2,s1,s2,m1,m2)
    s = sqrt(((n1-1)*s1^2+(n2-1)*s2^2)/(n1+n2-1)+n1*n2*(m1-m2)^2/(n1+n2)/(n1+n2-1));
end

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
