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


