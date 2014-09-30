%% HW1 Problem 2: Keplerian Orbital Elements to  Cartesian Coordinates
fprintf('\n');
clearvars -except function_list pub_opt
close all

r = [-2436.45; -2436.45; 6891.037]; % km
v = [5.088611; -5.088611; 0.0]; % km/s
state = [r;v];
oe = cart2oe(state);
new_state = oe2cart(oe);
state_diff = new_state - state;

fprintf('Computed State Vector\n')
new_state
fprintf('Delta between original state vector and computed vector')
state_diff