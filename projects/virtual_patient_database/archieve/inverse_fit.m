clear; close all;

file_name = "output.txt";
mmHg_to_Pa = 133.322368;
p0 = 1e5;
m3s_to_mls = 1e6;

raw_data = importdata(file_name);

data = raw_data.data;
cases = raw_data.textdata;

p_rad_d = (data(:,1)-p0)/mmHg_to_Pa;
p_rad_s = (data(:,2)-p0)/mmHg_to_Pa;
p_aor_d = (data(:,3)-p0)/mmHg_to_Pa;
p_aor_s = (data(:,4)-p0)/mmHg_to_Pa;
p_car_d = (data(:,5)-p0)/mmHg_to_Pa;
p_car_s = (data(:,6)-p0)/mmHg_to_Pa;
q_fem_a_l = data(:,7)*m3s_to_mls;
q_fem_a_r = data(:,8)*m3s_to_mls;

% values from literature
p_rad_d_fit = [65,12.84];
p_rad_s_fit = [123,22.75];
p_aor_d_fit = [65,7.14];
p_aor_s_fit = [103,8.24];
p_car_d_fit = [75.58,6.01];
p_car_s_fit = [122.78,12.01];
q_fem_a_fit = [5.84,2.11];


figure('position',[400,100,800,800]);
subplot(4,2,1);
x = (p_rad_d_fit(1)-3*p_rad_d_fit(2)):p_rad_d_fit(2)*.01:(p_rad_d_fit(1)+3*p_rad_d_fit(2));
y = normpdf(x,p_rad_d_fit(1),p_rad_d_fit(2));
plot(x,y);
grid on;

subplot(4,2,2);
x = (p_rad_s_fit(1)-3*p_rad_s_fit(2)):p_rad_s_fit(2)*.01:(p_rad_s_fit(1)+3*p_rad_s_fit(2));
y = normpdf(x,p_rad_s_fit(1),p_rad_s_fit(2));
plot(x,y);
grid on;

subplot(4,2,3);

subplot(4,2,4);

subplot(4,2,5);

subplot(4,2,6);

subplot(4,2,7);

subplot(4,2,8);


