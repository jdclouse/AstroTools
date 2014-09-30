%% HW1 Problem 1: Cartesian Coordinates to Keplerian Orbital Elements
fprintf('\n');
clearvars -except function_list pub_opt
close all

r = [-2436.45; -2436.45; 6891.037]; % km
v = [5.088611; -5.088611; 0.0]; % km/s
state = [r;v];
oe = cart2oe(state);
fprintf('a = %f km\n', oe(1))
fprintf('e = %f\n', oe(2))
fprintf('i = %f degrees\n', oe(3)*180/pi)
fprintf('RAAN = %f degrees\n', oe(4)*180/pi)
fprintf('Arg of Periapse = %f degrees\n', oe(5)*180/pi)
fprintf('True Anomaly = %f degrees\n', oe(6)*180/pi)