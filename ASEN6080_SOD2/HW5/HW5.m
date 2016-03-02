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
tic
storage = SRIF(state, P, meas_store, filter_opts);
toc

plot_cov_err_envelope(storage.cov_store(1:3,:), storage.state_store - true_state*1e3)
residual_plot(storage.pfr_store, sigs,'Q=0')


%% CKF

filter_opts.use_EKF = 0;
filter_opts.use_SNC = 0;

CKF_storage = KalmanFilter(state, meas_store, filter_opts);
plot_cov_err_envelope(CKF_storage.cov_store(1:3,:), CKF_storage.state_store - true_state*1e3)
residual_plot(CKF_storage.pfr_store, sigs,'Q=0')