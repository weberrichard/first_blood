clear;

l = 0.1;
n = 21;
x = linspace(0,l,n);
dx = x(2)-x(1);

v0 = 2;
v = v0*ones(n,1);
fi = zeros(n,1);
o2 = zeros(n,1);
fi_old = zeros(n,1);

fi_c_fut = []; % capillaries
fi_c = [];
o2_c = [];
dt_c = .1;
c_c = .9;

fi_p_fut = []; % pulmonary
fi_p = [];
o2_p = [];
dt_p = .1;
c_p = 1.;

dt_h = .1;

CFL = 1.;
dt = CFL*dx/v0;

fi0 = 5e12; % initial RBC count
fin = .2;
t_max = 5;

HG_RBC = 270e6; % number of hemoglobin per red blood cell

figure;
t=0;
while(t(end)<t_max)
%     if(t>5*dt)
%         v = -v;
%     end
    
    for i=2:n-1
        if(v(i)>0)
            fi(i) = fi_old(i) - v(i)*dt/dx*(fi_old(i)-fi_old(i-1));
        else
            fi(i) = fi_old(i) - v(i)*dt/dx*(fi_old(i)-fi_old(i+1));
        end
    end

    % downstream BC
    if(v(n)>0)
        fi(n) = fi_old(n) - v(n)*dt/dx*(fi_old(n)-fi_old(n-1));
    else
        fi(n) = fin;
    end

    % capillaries
    fi_c_fut = [fi_c_fut,fi(n)];
    if(t(end)>dt_c)
        fi_c = [fi_c,fi_c_fut(end-dt_c/dt)];
    else
        fi_c = [fi_c,fi0];
    end

    % pulmonary
    fi_p_fut = [fi_p_fut,fi_c];
    if(t(end)>dt_c+dt_p)
        fi_p = [fi_p,fi_p_fut(end-dt_p/dt)];
    else
        fi_p = [fi_p,fi0];
    end

    % upstream BC
    if(t(end)>dt_c+dt_p+dt_h)
        fi(1) = fi_p(end-dt_h/dt);
    else
        fi(1) = fi0;
    end
    
    % plotting stuff
    subplot(3,1,1);cla;
    plot(x,fi);
    hold on; grid on;
    plot(x,fi_old,'--');
    title(["t= ",num2str(t(end))]);
    subplot(3,1,2);
    plot(t,fi_c);
    subplot(3,1,3);
    plot(t,fi_p);

    pause(.1);

    fi_old = fi;
    t = [t,t(end)+dt];
end
