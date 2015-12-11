%% HW4 Problem 3
fprintf('\n');
clearvars -except function_list hw_pub toolsPath 
close all
CelestialConstants; % import useful constants

lat = 40.01*pi/180; % rad
lon = 254.83*pi/180; % rad
h = 1.615; %km
theta_GST = 103*pi/180; % rad

r_eci = ecef2eci(latlonalt2ECEF(lat, lon, h), theta_GST);
fprintf('r_eci = %.2f \n', r_eci(1));
fprintf('        %.2f km\n', r_eci(2));
fprintf('        %.2f \n', r_eci(3));