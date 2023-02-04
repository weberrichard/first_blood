clear;

% from McEniery2005
m = [172 178 183 258 429 430 280 39];
f = [133 101 165 301 495 509 290 38];
n = m+f;
sys_m_m = [103 105 109 113 115 117 118 120];
sys_m_std = [8 8 9 9 9 9 9 8];
sys_f_m = [98 101 105 109 115 118 119 120];
sys_f_std = [9 9 11 11 11 10 9 11];
dia_m_m = [74 75 78 79 80 78 76 75];
dia_m_std = [9.4 10.0 10.8 10.8 11.4 11.4 11.4 12.0];
dia_f_m = [73 74 75 76 77 75 63 71];
dia_f_std = [10.8 11.4 13.6 13.6 13.6 12.8 12.0 16.3];

sys_m_mean = sum(m.*sys_m_m)/sum(m);
sys_f_mean = sum(f.*sys_f_m)/sum(f);

dia_m_mean = sum(m.*dia_m_m)/sum(m);
dia_f_mean = sum(f.*dia_f_m)/sum(f);

sys_m_std = sdc_gen(sys_m_m,m,sys_m_std);
sys_f_std = sdc_gen(sys_f_m,f,sys_f_std);
dia_m_std = sdc_gen(dia_m_m,m,dia_m_std);
dia_f_std = sdc_gen(dia_f_m,f,dia_f_std);

sys_m = (sum(m)*sys_m_mean+sum(f)*sys_f_mean)/(sum(m)+sum(f));
dia_m = (sum(m)*dia_m_mean+sum(f)*dia_f_mean)/(sum(m)+sum(f));

sys_std = sdc_gen([sys_m_mean,sys_f_mean],[sum(m),sum(f)],[sys_m_std,sys_f_std]);
dia_std = sdc_gen([dia_m_mean,dia_f_mean],[sum(m),sum(f)],[dia_m_std,dia_f_std]);

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