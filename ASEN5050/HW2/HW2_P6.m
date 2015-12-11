%% HW2 Problem 6
fprintf('\n');
clearvars -except function_list hw_pub toolsPath 
close all
CelestialConstants; % import useful constants

a = 0.387; %AU
e = 0.205;
TU2days = 365.25/(2*pi);

P = 2*pi*sqrt(a*a*a)*TU2days; %s
ra = a*(1+e);
rp = a*(1-e);

vp = sqrt(2/rp - 1/a)*au2km/TU2days/day2sec;
fprintf('Aphelion = %.3f AU\n',ra);
fprintf('Perihelion = %.3f AU\n',rp);
fprintf('Perihelion speed = %.1f km/s\n',vp);