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
filter_opts.propagator_opts.J3.use = 0;
storage = UnscentedKalmanFilter(state, P, meas_store, filter_opts);

% pos_cov = arrayfun(@(ii) norm(storage.Pt_store(:,ii)), 1:num_obs);
% plot_cov_err_envelope(pos_cov, storage.X_est_store - true_state*1e3)
plot_cov_err_envelope(storage.Pt_store(1:3,1:500), storage.X_est_store - true_state(:,1:500)*1e3)
% plot_cov_err_envelope(storage.Pt_store(1:3,:), storage.X_est_store - true_state*1e3)
title('UKF State Error, with covariance envelope')

residual_plot(storage.pfr_store, sigs)
title('UKF, No Q')
storage.range_RMS = sqrt(sum(storage.pfr_store(1,:).*storage.pfr_store(1,:))/num_obs);
storage.rangerate_RMS = sqrt(sum(storage.pfr_store(2,:).*storage.pfr_store(2,:))/num_obs);