%% HW5

%% Initialize
% clearvars -except function_list 
global function_list;
function_list = {};
addpath('C:\Users\John\Documents\Astro\ASEN6080_SOD2\tools\')

% Load the params file
filter_params;

%% 
num_obs = length(meas_store);
sig_range = 0.01; % m
sig_rangerate = 0.001; %m/s
sigs = [sig_range; sig_rangerate];

% Using UKF
% P = eye(6)*1e2;
% filter_opts.use_SNC = 0;
filter_opts.propagator_opts.J3.use = 0;
filter_opts.ode_opts = odeset('RelTol', 1e-12, 'AbsTol', 1e-20);
[unique_meas_times, i_times, i_uniques] = unique(meas_store(:,2),'stable');
[~,X] = ode45(@two_body_state_dot, unique_meas_times, [state; reshape(eye(6),36,1)], ...
            filter_opts.ode_opts, filter_opts.propagator_opts);
filter_opts.integrate_ref_state = 0;
filter_opts.ref_state = X(i_uniques,:);
filter_opts.ref_state = X;
%% 
tic
storage = SRIF(state, P, meas_store, filter_opts);
toc

% plot_cov_err_envelope(storage.cov_store(1:3,:), storage.state_store - true_state*1e3)
residual_plot(storage.pfr_store, sigs,'SRIF, R not Triangular')


%% CKF
% 
% filter_opts.use_EKF = 0;
% filter_opts.use_SNC = 0;
% filter_opts.use_DMC = 0;
% 
% CKF_storage = KalmanFilter(state, meas_store, filter_opts);
% plot_cov_err_envelope(CKF_storage.cov_store(1:3,:), CKF_storage.state_store - true_state*1e3)
% residual_plot(CKF_storage.pfr_store, sigs,'CKF')