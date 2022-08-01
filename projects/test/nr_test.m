clear;

n = 4;
err = 1e-4;

x = ones(n,1);
dx = 1e10;
d = 1e-2;

while(norm(dx)>err)
    f0 = test3(x);
    Jac = [[4*x(1),3,4,-5];[1,-1,1,-1];[-3,4,2,1];[1,0,1,0]];
%     Jac = 6*x;
    
    dx = Jac\-f0;
    x = x+dx;
    norm(dx)
end


function out = test(x)
    A = [[2,3,4,-5];[1,-1,1,-1];[-3,4,2,1];[1,0,3*x(3),0]];
    b = [2,0,-4,1]';
    
    out = A*x-b;
end


function out = test2(x)
    A = 3*x;
    b = 1;
    
    out = A*x-b;
end

function out = test3(x)
    out(1,1) = 2*x(1)^2+3*x(2)+4*x(3)-5*x(4)-2;
    out(2,1) = 1*x(1)-1*x(2)+1*x(3)-1*x(4)-0;
    out(3,1) =-3*x(1)+4*x(2)+2*x(3)+1*x(4)+4;
    out(4,1) = 1*x(1)+0*x(2)+1*x(3)-0*x(4)-1;
end