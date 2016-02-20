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
sig_range = 0.01; % m
sig_rangerate = 0.001; %m/s
sigs = [sig_range; sig_rangerate];

% Using UKF
P = eye(6)*1e2;
filter_opts.use_SNC = 0;
storage = UnscentedKalmanFilter(state, P, meas_store, filter_opts);

pos_cov = arrayfun(@(ii) norm(storage.Pt_store(:,ii)), 1:num_obs);
plot_cov_err_envelope(pos_cov, storage.X_est_store - true_state*1e3)
title('UKF State Error, with covariance envelope')

residual_plot(storage.pfr_store, sigs)
title('UKF, No Q')

%% Adding process noise
filter_opts.use_SNC = 1;
filter_opts.SNC_Q = eye(3)*1e-5;
storage_withQ = UnscentedKalmanFilter(state, P, meas_store, filter_opts);

pos_cov = arrayfun(@(ii) norm(storage_withQ.Pt_store(:,ii)), 1:num_obs);
plot_cov_err_envelope(pos_cov, storage_withQ.X_est_store - true_state*1e3)
title('UKF State Error, with covariance envelope')

residual_plot(storage_withQ.pfr_store, sigs)
title('UKF, Q = 1e-5')

filter_opts.SNC_Q = eye(3)*1e-6;
storage_withQ2 = UnscentedKalmanFilter(state, P, meas_store, filter_opts);

pos_cov = arrayfun(@(ii) norm(storage_withQ2.Pt_store(:,ii)), 1:num_obs);
plot_cov_err_envelope(pos_cov, storage_withQ2.X_est_store - true_state*1e3)
title('UKF State Error, with covariance envelope')

residual_plot(storage_withQ2.pfr_store, sigs)
title('UKF, Q = 1e-6')
