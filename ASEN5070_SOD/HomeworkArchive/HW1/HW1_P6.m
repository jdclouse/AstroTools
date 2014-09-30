%% HW1 Problem 6: Numerical estimation of initial state
fprintf('\n');
clearvars -except function_list pub_opt
close all

state = [1.5
    10
    2.2
    0.5
    0.3];

station_loc = [1.0; 1.0;];

obs = [7.0, 0
    8.00390597, 1
    8.94427191, 2
    9.801147892, 3
    10.630145813, 4];
    
tol = 1e-6;

delta = 1;
while delta > tol
    jac = compute_cost_fcn_jacobian(obs, state, station_loc);
    cost_fnc = compute_obs_2d(obs, state, station_loc);
    state_est = state - jac\cost_fnc;
    delta = norm(state_est - state);
    state = state_est;
end

fprintf('True Initial State:\n');
state