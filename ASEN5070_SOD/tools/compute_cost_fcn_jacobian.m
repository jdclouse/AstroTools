function jac = compute_cost_fcn_jacobian(obs, state, station_loc)
%compute_cost_fcn_jacobian   Compute the cost function jacobian for a 2D
%flat-earth problem.
fcnPrintQueue(mfilename('fullpath')) % Add this code to code appendix

num_unks = length(state);
jac = zeros(num_unks, num_unks);
for ii = 1:num_unks
    jac(ii, :) = compute_range_partials(obs(ii,:), state, station_loc);
end

end

function row = compute_range_partials(observation, state, station_loc)
%compute_range_partials   Compute a row of the flat-earth-problem jacobian.
X0 = state(1);
Y0 = state(2);
X0_dot = state(3);
Y0_dot = state(4);
g = state(5);
Xs = station_loc(1);
Ys = station_loc(2);

rho = observation(1);
t = observation(2);

xterm = -(X0 - Xs + X0_dot*t)/rho;
yterm = -(Y0 - Ys + Y0_dot*t - g*t*t/2)/rho;

row = [xterm, yterm, xterm*t, yterm*t, -yterm*t*t/2];

end