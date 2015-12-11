%% HW4 Problem 2
fprintf('\n');
clearvars -except function_list hw_pub toolsPath 
close all
CelestialConstants; % import useful constants

r_eci = [-5634 -2645 2834]'; % km
theta_GST = 82.75*pi/180;

[lat, lon, alt] = ECEF2latlonalt(eci2ecef(r_eci, theta_GST));
fprintf('Latitude: %.2f deg\n', lat*180/pi);
fprintf('Longitude: %.2f deg\n', lon*180/pi);
fprintf('Altitude: %.2f km\n', alt);

% ecef2eci(latlonalt2ECEF(lat, lon, alt), theta_GST)