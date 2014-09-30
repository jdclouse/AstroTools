function cost_fnc = compute_obs_2d(obs, state, station_loc)
%compute_obs_2d   Compute the cost function given the est. state and
%   observations.  Currently 2D, range observations, num measurements are
%   exactly the num unks
fcnPrintQueue(mfilename('fullpath')) % Add this code to code appendix

num_unks = length(state);
cost_fnc = zeros(num_unks, 1);
for ii = 1:num_unks
    cost_fnc(ii) = obs(ii,1) ...
        - compute_range(state, station_loc, obs(ii,2));
end

end

function rho = compute_range(state, station_loc, t)
%compute_range   Compute the calculated range, given state.
X0 = state(1);
Y0 = state(2);
X0_dot = state(3);
Y0_dot = state(4);
g = state(5);
Xs = station_loc(1);
Ys = station_loc(2);

xterm = (X0-Xs+X0_dot*t)*(X0-Xs+X0_dot*t);
yterm = (Y0-Ys+Y0_dot*t-g*t*t/2)*(Y0-Ys+Y0_dot*t-g*t*t/2);
rho = sqrt(xterm + yterm);

end

