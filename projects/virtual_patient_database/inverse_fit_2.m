clear; close all;

file_name = "output.txt";
mmHg_to_Pa = 133.322368;
p0 = 1e5;
m3s_to_mls = 1e6;

raw_data = importdata(file_name);

meas = raw_data.data;
cases = raw_data.textdata;

% converting simulation data
meas(:,1:6) = (meas(:,1:6)-p0)/mmHg_to_Pa;
meas(:,7:8) = meas(:,7:8)*m3s_to_mls;

% values from literature
lite= [65,12.84; % radial dias pres mu, std
123,22.75; % radial sys pres mu, std
65,7.14; % aorta dias pres mu, std
103,8.24; % aorta sys pres mu, std
75.58,6.01; % carotis dias pres mu, std
122.78,12.01; % carotis sys pres mu, std
5.84,2.11; % femoral volume flow rate average mu, std
5.84,2.11]; % femoral volume flow rate average mu, std
    
% plotting stuff
figure('position',[400,100,800,800]);

% titles
titles = ["Radial diastole","Radial systole","Aorta diastole","Aorta systole","Carotis diastole","Carotis systole","Femoral right","Femoral left"];

% plotting literature data
for i=1:8
   subplot(4,2,i);
   x = lite(i,1)-3*lite(i,2):lite(i,2)*.01:lite(i,1)+3*lite(i,2);
   y = normpdf(x,lite(i,1),lite(i,2));
   plot(x,y,'linewidth',1.5);   
   grid on; hold on;
   if(i<7)
      xlabel("pressure [mmHg]"); 
   else
      xlabel("volume flow rate [ml/s]");
   end
   ylabel("probability");
   title(titles(i));
end

% plotting the original simulation results
for i=1:8
    subplot(4,2,i);
    histogram(meas(:,i),'Normalization','pdf');
end

% original fitness function
n = size(meas,1);
m = size(meas,2);

xhi2orig = zeros(m,1);
for i=1:m
    r = 2*round(log(m)/log(2)+1);
    p = prctile(meas(:,i),linspace(0,100,r));
    p(1) = p(1) - 0.001;
    nu = sum(meas(:,i)>p(1:end-1) & meas(:,i)<=p(2:end),1)/(n-1);
    p(1) = p(1) + 0.001;
    pi = normcdf(p,lite(i,1),lite(i,2));
    pi(1) = 0;
    pi(end) = 1;
    pi = pi(2:end)-pi(1:end-1);
    xhi2orig(i) = sum((pi-nu).^2.);
end
    
% selecting the database
N = 1;
meas_fit = meas;
for k = 1:N
    n = size(meas_fit,1);
    good_index = zeros(n,1);
    for j=1:n
        % fitness function
        xhi2 = zeros(m,1);
        for i=1:m
            r = 2*round(log(m)/log(2)+1);
            p = prctile(meas_fit([1:j-1,j+1:end],i),linspace(0,100,r));
            p(1) = p(1) - 0.001;
            nu = sum(meas_fit([1:j-1,j+1:end],i)>p(1:end-1) & meas_fit([1:j-1,j+1:end],i)<=p(2:end),1)/(n-1);
            p(1) = p(1) + 0.001;
            pi = normcdf(p,lite(i,1),lite(i,2));
            pi(1) = 0;
            pi(end) = 1;
            pi = pi(2:end)-pi(1:end-1);
            xhi2(i) = sum((pi-nu).^2.);
        end
%         if(sum(xhi2)>sum(xhi2orig))
        if(mean(xhi2)>mean(xhi2orig))
           good_index(j) = 1; 
        end
    end
    meas_fit = meas(logical(good_index),:);

    % plotting the original simulation results
    for i=1:8
        subplot(4,2,i);
        histogram(meas_fit(:,i),'Normalization','pdf');
        if(i>=7)
           xlim([-20,50]); 
        end
    end
end