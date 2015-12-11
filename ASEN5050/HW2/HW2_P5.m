%% HW2 Problem 5
fprintf('\n');
clearvars -except function_list hw_pub toolsPath 
close all
CelestialConstants; % import useful constants

r = 8200; %km
vcirc = sqrt(Earth.mu/r);

fprintf('Satellite speed should be %.3f km/s\n', vcirc)
fprintf('Satellite FPA should be 0 degrees\n')