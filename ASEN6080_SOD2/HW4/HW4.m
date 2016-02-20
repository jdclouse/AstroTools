%% Filter Data

%% Initialize
% clearvars -except function_list 
global function_list;
function_list = {};
addpath('C:\Users\John\Documents\Astro\ASEN5070_SOD\tools\')
addpath('C:\Users\John\Documents\Astro\ASEN5050\tools\')

% Load the params file
filter_params;

%% 
sig_store = [];
rms_store = [];
pos_RMS = [];
vel_RMS = [];
% skip_obs = 100; % observations to skip when computing error RMS.
num_obs = length(meas_store);

% Using UKF
P = eye(6)*1e2;
storage = UnscentedKalmanFilter(state, P, meas_store, filter_opts);

plot_cov_err_envelope(storage.cov_store, storage.state_store - true_state*1e3)
title('CKF State Error, with covariance envelope')
state_error = storage.smoothed_state_store - true_state*1e3;
plot_cov_err_envelope(storage.P_smoothed_diag, state_error)
title('CKF Smoothed State Error, with smoothed covariance envelope')

