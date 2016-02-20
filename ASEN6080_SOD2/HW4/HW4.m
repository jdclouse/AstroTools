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
filter_opts.use_SNC = 0;
storage = UnscentedKalmanFilter(state, P, meas_store, filter_opts);

pos_cov = arrayfun(@(ii) norm(storage.Pt_store(:,ii)), 1:num_obs);
plot_cov_err_envelope(pos_cov, storage.X_est_store - true_state*1e3)
title('UKF State Error, with covariance envelope')

%% Adding process noise
filter_opts.use_SNC = 1;
storage_withQ = UnscentedKalmanFilter(state, P, meas_store, filter_opts);

pos_cov = arrayfun(@(ii) norm(storage_withQ.Pt_store(:,ii)), 1:num_obs);
plot_cov_err_envelope(pos_cov, storage_withQ.X_est_store - true_state*1e3)
title('UKF State Error, with covariance envelope')

sig_range = 0.01; % m
sig_rangerate = 0.001; %m/s
figure
subplot(2,1,1)
plot(1:num_obs, storage.pfr_store(1,:),'.','LineWidth',1)
hold on
plot(1:num_obs,3*sig_range*ones(1,num_obs),'r--')
plot(1:num_obs,-3*sig_range*ones(1,num_obs),'r--')
% title(sprintf('Range RMS = %.4e m',output.range_RMS))
ylabel('m')
subplot(2,1,2)
plot(1:num_obs, storage.pfr_store(2,:),'.','LineWidth',1)
hold on
plot(1:num_obs,3*sig_rangerate*ones(1,num_obs),'r--')
plot(1:num_obs,-3*sig_rangerate*ones(1,num_obs),'r--')
% title(sprintf('Range-Rate RMS = %.4e m/s',output.rangerate_RMS))
ylabel('m/s'),xlabel('Observation')