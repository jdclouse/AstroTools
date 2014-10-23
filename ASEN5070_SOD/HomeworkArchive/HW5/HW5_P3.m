%% HW5 Problem 3
%% Initialize
fprintf('\n');
clearvars -except function_list pub_opt
close all
% Bring in answers to compare
hw5_p1_answers

%% Estimate the initial state
rho = [6.37687486186586
    5.50318198665912
    5.94513302809067
    6.30210798411686
    5.19084347133671
    6.31368240334678
    5.80399842220377
    5.45115048359871
    5.91089305965839
    5.67697312013520
    5.25263404969825];
rho_dot = [-0.00317546143535849
     1.17587430814596
    -1.47058865193489
     0.489030779000695
     0.993054430595876
    -1.40470245576321
     0.939807575607138
     0.425908088320457
    -1.47604467619908
     1.42173765213734
    -0.12082311844776];

Y = [rho'; rho_dot'];
W = inv([0.0625 0; 0 0.01]);
x=[0; 0];
W_bar = inv([1000 0; 0 100]);
x_bar = x;

[state_est, BLS_info] = BLS_spring( x, Y, W, x_bar, W_bar );
fprintf('x0 is %.4f m\n', state_est(1))
fprintf('v0 is %.4f m/s\n', state_est(2))
fprintf('Range RMS is %.3f m\n', BLS_info.RMS(1))
fprintf('Range rate RMS is %.4f m/s\n', BLS_info.RMS(2))
sig_x = sqrt(BLS_info.P0(1,1));
sig_v = sqrt(BLS_info.P0(2,2));
rho_xv=BLS_info.P0(1,2)/sig_x/sig_v;
fprintf('x standard deviation is %.4f m/s\n', sig_x)
fprintf('v standard deviation is %.4f m/s\n', sig_v)
fprintf('correlation of x and v is %.4f m/s\n', rho_xv)