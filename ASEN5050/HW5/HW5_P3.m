%% HW5 Problem 3
fprintf('\n');
clearvars -except function_list hw_pub toolsPath 
close all
CelestialConstants; % import useful constants

phase_angle = 30*pi/180;
h = 6e3; %km
a = Earth.R + h;
n = sqrt(Earth.mu/a^3);
t_phase = (2*pi + phase_angle)/n
t_phase/3600;
a_phase = ((t_phase/(2*pi))^2*Earth.mu)^(1/3);
dv1 = sqrt(2*Earth.mu/a-Earth.mu/a_phase) - sqrt(Earth.mu/a)
dv2 = sqrt(Earth.mu/a) - sqrt(2*Earth.mu/a-Earth.mu/a_phase)
dv_tot = abs(dv1) + abs(dv2)