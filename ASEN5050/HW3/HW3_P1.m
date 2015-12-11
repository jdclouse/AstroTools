%% HW3 Problem 1
fprintf('\n');
clearvars -except function_list hw_pub toolsPath 
close all
CelestialConstants; % import useful constants

r = [-5633.9; -2644.9; 2834.4];
v = [2.425; -7.103; -1.800];

[a,e,i,RAAN,w,f] = cart2OE(r,v,Earth.mu);
fprintf('Semi-major axis: %.2f km\n',a);
fprintf('Eccentricity: %.5f \n',e);
fprintf('Inclination: %.2f deg\n',i * 180/pi);
fprintf('Right Ascention of Asc. Node: %.2f deg\n',RAAN * 180/pi);
fprintf('Argument of Periapse: %.2f deg\n',w * 180/pi);
fprintf('True anomaly: %.2f deg\n',f * 180/pi);