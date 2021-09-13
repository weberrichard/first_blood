clear;
% Wilson's example
b=[32,23,33,31]';
A=[10,7,8,7;7,5,6,5;8,6,10,9;7,5,9,10];

C = cond(A);

% P = eye(4);
% P(1,1) = 10;
[P,~] = eig(A);

A2 = (inv(P)*A*P);

C2 = cond(A2);

A3 = A/norm(A);

