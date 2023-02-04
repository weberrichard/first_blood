clear;

qm = 44; % ml/min/m2
qs = 6;

Dm = 0.436; % cm
Ds = 0.013;

% conv to SI
% qm = qm/1000/60;
% qs = qm/1000/60;

Dm = Dm/100;
Ds = Ds/100;

Qm = qm*Dm^2*pi/4;

vm = 9.8; % cm/s
vm = vm/100;
Qm2 = vm*Dm^2*pi/4;
Qm2 = Qm2*1000*60;

% conv to ml/min
% Qm = Qm*1000*60;