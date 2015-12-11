%% HW3 Problem 2
fprintf('\n');
clearvars -except function_list hw_pub toolsPath 
close all
CelestialConstants; % import useful constants

p = 11067.790;
e=0.83285;
a = p/(1-e*e);
i = 87.87*pi/180;
RAAN = 227.89*pi/180;
w = 53.38*pi/180;
f = 92.335*pi/180;
[r, v ] = OE2cart( a,e,i,RAAN,w,f,Earth.mu);

fprintf('r_ECI = %.2f\n',r(1));
fprintf('        %.2f km\n',r(2));
fprintf('        %.2f\n\n',r(3));
fprintf('v_ECI = %.3f\n',v(1));
fprintf('        %.3f km/s\n',v(2));
fprintf('        %.3f\n',v(3));