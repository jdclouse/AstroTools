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
storage.range_RMS = sqrt(sum(storage.pfr_store(1,:).*storage.pfr_store(1,:))/num_obs);
storage.rangerate_RMS = sqrt(sum(storage.pfr_store(2,:).*storage.pfr_store(2,:))/num_obs);

%% Adding process noise
filter_opts.use_SNC = 1;
filter_opts.SNC_Q = eye(3)*1e-5;
storage_withQ = UnscentedKalmanFilter(state, P, meas_store, filter_opts);

pos_cov = arrayfun(@(ii) norm(storage_withQ.Pt_store(:,ii)), 1:num_obs);
plot_cov_err_envelope(pos_cov, storage_withQ.X_est_store - true_state*1e3)
title('UKF State Error, with covariance envelope')

residual_plot(storage_withQ.pfr_store, sigs)
title('UKF, Q = 1e-5')
storage_withQ.range_RMS = sqrt(sum(storage_withQ.pfr_store(1,:).*storage_withQ.pfr_store(1,:))/num_obs);
storage_withQ.rangerate_RMS = sqrt(sum(storage_withQ.pfr_store(2,:).*storage_withQ.pfr_store(2,:))/num_obs);

filter_opts.SNC_Q = eye(3)*1e-6;
storage_withQ2 = UnscentedKalmanFilter(state, P, meas_store, filter_opts);

pos_cov = arrayfun(@(ii) norm(storage_withQ2.Pt_store(:,ii)), 1:num_obs);
plot_cov_err_envelope(pos_cov, storage_withQ2.X_est_store - true_state*1e3)
title('UKF State Error, with covariance envelope')

residual_plot(storage_withQ2.pfr_store, sigs)
title('UKF, Q = 1e-6')
storage_withQ2.range_RMS = sqrt(sum(storage_withQ2.pfr_store(1,:).*storage_withQ2.pfr_store(1,:))/num_obs);
storage_withQ2.rangerate_RMS = sqrt(sum(storage_withQ2.pfr_store(2,:).*storage_withQ2.pfr_store(2,:))/num_obs);

%% Reduce alpha
filter_opts.UKF.alpha = 1e-4; 
storage_reduced_alpha = UnscentedKalmanFilter(state, P, meas_store, filter_opts);

pos_cov = arrayfun(@(ii) norm(storage_reduced_alpha.Pt_store(:,ii)), 1:num_obs);
plot_cov_err_envelope(pos_cov, storage_reduced_alpha.X_est_store - true_state*1e3)
title('UKF State Error, with covariance envelope')

residual_plot(storage_reduced_alpha.pfr_store, sigs)
title('UKF, alpha = 1e-4')
storage_reduced_alpha.range_RMS = sqrt(sum(storage_reduced_alpha.pfr_store(1,:).*storage_reduced_alpha.pfr_store(1,:))/num_obs);
storage_reduced_alpha.rangerate_RMS = sqrt(sum(storage_reduced_alpha.pfr_store(2,:).*storage_reduced_alpha.pfr_store(2,:))/num_obs);

%% Large a priori error
filter_opts.UKF.alpha = 1; 
state_big_ap_error = state_i*1e3 + [500;500;500;10;10;10];
P = eye(6)*500;
storage_big_ap = UnscentedKalmanFilter(state_big_ap_error, P, meas_store, filter_opts);

pos_cov = arrayfun(@(ii) norm(storage_big_ap.Pt_store(:,ii)), 1:num_obs);
plot_cov_err_envelope(pos_cov, storage_big_ap.X_est_store - true_state*1e3)
title('UKF State Error, with covariance envelope')

residual_plot(storage_big_ap.pfr_store, sigs)
title('UKF, large a priori error')
storage_big_ap.range_RMS = sqrt(sum(storage_big_ap.pfr_store(1,:).*storage_big_ap.pfr_store(1,:))/num_obs);
storage_big_ap.rangerate_RMS = sqrt(sum(storage_big_ap.pfr_store(2,:).*storage_big_ap.pfr_store(2,:))/num_obs);


storage_EKFbig_ap = KalmanFilter(state_big_ap_error, meas_store, filter_opts);
plot_cov_err_envelope(pos_cov, storage_EKFbig_ap.state_store - true_state*1e3)
title('EKF State Error, with covariance envelope')

residual_plot(storage_EKFbig_ap.pfr_store, sigs)
title('EKF, large a priori error')

%% Add J3
P = eye(6)*1e2;
filter_opts.propagator_opts.J3.use = 1;
storage_J3 = UnscentedKalmanFilter(state, P, meas_store, filter_opts);
pos_cov = arrayfun(@(ii) norm(storage_J3.Pt_store(:,ii)), 1:num_obs);
plot_cov_err_envelope(pos_cov, storage_J3.X_est_store - true_state*1e3)
title('UKF State Error, with covariance envelope')

residual_plot(storage_J3.pfr_store, sigs)
title('UKF, with J3')