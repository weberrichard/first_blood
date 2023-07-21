clear;

l = 0.1;
n = 21;
x = linspace(0,l,n);
dx = x(2)-x(1);

v0 = 1;
v = v0*ones(n,1);
fi = zeros(n,1);
fi_old = zeros(n,1);

CFL = 1.;
dt = CFL*dx/v0;
fi0 = .5;
fin = .2;
t_max = 5;

figure;
t=0;
while(t<t_max)
    if(t>5*dt)
        v = -v;
    end
    
    for i=2:n-1
        if(v(i)>0)
            fi(i) = fi_old(i) - v(i)*dt/dx*(fi_old(i)-fi_old(i-1));
        else
            fi(i) = fi_old(i) - v(i)*dt/dx*(fi_old(i)-fi_old(i+1));
        end
    end
    
    % upstream BC
    if(v(1)>0)
        fi(1) = fi0;
    else
        fi(1) = fi_old(1) - v(1)*dt/dx*(fi_old(1)-fi_old(2));
    end

    % downstream BC
    if(v(n)>0)
        fi(n) = fi_old(n) - v(n)*dt/dx*(fi_old(n)-fi_old(n-1));
    else
        fi(n) = fin;
    end
    
    % plotting stuff
    cla;
    plot(x,fi);
    hold on; grid on;
    plot(x,fi_old,'--');
    
    fi_old = fi;
end
