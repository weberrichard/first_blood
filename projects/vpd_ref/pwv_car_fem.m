clear;

% aorta -> carotis
a1 = 7;
L1 = .094;

% aorta -> femoral
a2 = 8.5;
L2 = 0.8;

dt01 = L1/a1;
dt02 = L2/a2;
dt12 = dt02-dt01;

a_orig = (L1+L2)/dt12;

LL = a2/a1*L1;

a_real = (L2-LL)/dt12;

a_cor = (L2-L1)/dt12;