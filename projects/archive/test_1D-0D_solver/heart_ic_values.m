clear; close all;

C_artial = 3.3003e-8;
C_system = 9.9758e-9;
C_elastance =  8.0016e-06;

p_artial = 1.0665e+05; %100842.3284;
p_ventri = 100407.96;
p_aorta =  1e5; %112000;

p0 = 1e5;

E0 = 0.5*C_artial*(p_artial-p0)^2 + 0.5*C_elastance*(p_ventri-p0)^2 + 0.5*C_system*(p_aorta-p0)^2;

p_artial_new = sqrt((2*E0-C_elastance*(p_ventri-p0)^2)/C_artial) + p0;
