%% HW5 Problem 2
fprintf('\n');
clearvars -except function_list hw_pub toolsPath 
close all
CelestialConstants; % import useful constants

a_xfer = (Earth.a + Mars.a)/2;
% Find the Lead Angle first:
n_mars = sqrt(Sun.mu/Mars.a^3);
T_xfer = pi*sqrt(a_xfer^3/Sun.mu);
lead_angle = n_mars*T_xfer;
phase_angle = lead_angle - pi;
phase_angle*180/pi

% Synodic period
n_earth = sqrt(Sun.mu/Earth.a^3);
syn_period = 2*pi/(n_earth-n_mars)/day2sec/365.25