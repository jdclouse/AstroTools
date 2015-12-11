%% HW2 Problem 3
fprintf('\n');
clearvars -except function_list hw_pub toolsPath 
close all
CelestialConstants; % import useful constants

hp_actual = 798; %km
ha_actual = 817; %km
hp_planned = 794; %km
ha_planned = 814; %km

a_actual = (hp_actual + ha_actual + 2*Earth.R)/2;
a_planned = (hp_planned + ha_planned + 2*Earth.R)/2;

e_actual = 1-(hp_actual+Earth.R)/a_actual;
e_planned = 1-(hp_planned+Earth.R)/a_planned;

P_actual = 2*pi*sqrt(a_actual*a_actual*a_actual/Earth.mu);
P_planned = 2*pi*sqrt(a_planned*a_planned*a_planned/Earth.mu);

fprintf('Semimajor axis:\n')
fprintf('\t%.2f km difference\n', abs(a_planned-a_actual))
fprintf('\t%.2f%% error\n\n', abs(a_planned-a_actual)/a_planned*100)
fprintf('Eccentricity:\n')
fprintf('\t%.6f difference\n', abs(e_planned-e_actual))
fprintf('\t%.2f%% error\n\n', abs(e_planned-e_actual)/e_planned*100)
fprintf('Period:\n')
fprintf('\t%.2f s difference\n', abs(P_planned-P_actual))
fprintf('\t%.2f%% error\n\n', abs(P_planned-P_actual)/P_planned*100)