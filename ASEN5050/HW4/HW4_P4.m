%% HW4 Problem 4
fprintf('\n');
clearvars -except function_list hw_pub toolsPath 
close all
CelestialConstants; % import useful constants

% Boulder, CO
lat = 40.01*pi/180; % rad
lon = 254.83*pi/180; % rad
h = 1.615; %km

r_sat_ecef = [-1681;-5173;4405];

topo = ecef2topo( r_sat_ecef, lat, lon, h);
fprintf('Azimuth: %.2f deg\n', topo(1)*180/pi);
fprintf('Elevation: %.2f deg\n', topo(2)*180/pi);
fprintf('Range: %.2f km\n', topo(3));