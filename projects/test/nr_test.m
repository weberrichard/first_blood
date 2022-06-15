clear;

n = 4;
err = 1e-4;

x = ones(n,1);
dx = 1e10;
d = 1e-2;

while(norm(dx)>err)
    f0 = test(x);
    Jac = [[2,3,4,-5];[1,-1,1,-1];[-3,4,2,1];[1,0,6*x(3),0]];
    
    dx = Jac\-f0;
    x = x+dx;
    norm(dx)
end


function out = test(x)

    A = [[2,3,4,-5];[1,-1,1,-1];[-3,4,2,1];[1,0,3*x(3),0]];
    b = [2,0,-4,1]';
    
    out = A*x-b;
end
