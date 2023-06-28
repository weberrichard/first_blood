
clear all , close all , clc
% set(0, ' defaultaxesfontsize ' ,12,' defaultlinelinewidth ' ,2)
L =100; D =0.4; lambda =0.02; a = 1300; rho =1000; vin =1; pout =1e5;
Npts =20; dx=L/( Npts -1); dt=dx/a; Tp =2*L/a; tend =100;
v= zeros(Npts ,1); p=1e5* ones(Npts ,1); vnew =v; pnew =p; alpha =v; beta =v; vends =[0 0 0];
% v_s = sqrt(( ptank - pout )/( lambda *L/D*rho /2) );
t=0;
while t< tend
for i=1: Npts
alpha (i)=p(i)+rho*a*v(i); beta(i) =p(i)-rho*a*v(i);
end
%% Update internal points with explicit scheme
for i=2: Npts -1
Sl=- lambda/2/D*v(i -1)* abs(v(i -1)); alpha_new = alpha(i -1)+dt* rho*a*Sl;
Sr=- lambda/2/D*v(i +1)* abs(v(i +1)); beta_new = beta(i +1) -dt* rho*a*Sr;
vnew(i)=( alpha_new - beta_new ) /2/ rho/a;
pnew(i)=( alpha_new + beta_new ) /2;
end
%% Left boundary : fixed pressure
vnew(1)= vin ; pnew(1) = rho*a*vnew(1) + beta(2) - dt*rho*a*lambda/2/D*v(1)*abs(v(1));
%% Right boundary : fixed pressure
pnew(end)= pout ; vnew(end)=(alpha(end-1)-dt*rho*a*lambda/2/D*v(end-1)*abs(v(end-1)) - pnew(end))/ rho/a;
%% Close step
t=t+dt; p= pnew ; v= vnew ; vends =[ vends ; t vnew(1) vnew(end)];
%% Plot distributions as movie
% figure(1) , xx= linspace(0,L, Npts );
% subplot(2 ,1 ,1) , plot(xx ,p/1e5)
% grid on , ylabel('p, bar ')
% title([ 't=',num2str(t),'s, ',num2str( round(t/ tend *100) ),'%'])
% subplot(2 ,1 ,2) , plot(xx ,v,'k' ,[0 L],vin *[1 1], 'r--')
% grid on , ylabel('v, m/s'), xlabel('x, m'), drawnow
end
%% Plot velocities @ x=0,L
figure(2)
tt= vends(: ,1);
plot(tt , vends(: ,2) ,'b',tt , vends(: ,3) ,'k' ,[0 t], vin *[1 1], 'r--')
xlabel('t, s'), ylabel ('v, m/s'), grid on
legend('v (0) ','v(L)','v steady ','Location ','best ')
