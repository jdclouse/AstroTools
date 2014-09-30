%% HW2 Problem 1: PRN generation
%% Initialize
fprintf('\n');
clearvars -except function_list pub_opt
close all

%% 
load('GPS2009L1L2_data.mat');
time = data(:,1);
lat = data(:,2);
long = data(:,3);
ell_h = data(:,4);
num_sat = data(:,5);
subplot(2,1,1)
plot(time, ell_h)
subplot(2,1,2)
plot(time, num_sat)

%% 2
mean_lat = mean(lat);
lat_err = lat - mean_lat;
mean_long = mean(long);
long_err = long - mean_long;

re = 6378e3; % m
NS_err = 2*pi*re*lat_err/360; 
EW_err = 2*pi*re*long_err/360; 

figure
plot(EW_err, NS_err, '.')
grid on
axis equal
hold on
nn = length(NS_err);
err_radii = sqrt((NS_err.*NS_err)+(EW_err.*EW_err));
CEP_50_r = err_radii(round(nn/2));
drawellipse(CEP_50_r,CEP_50_r,0,0,0,'r-.');

P = cov([EW_err NS_err]);

