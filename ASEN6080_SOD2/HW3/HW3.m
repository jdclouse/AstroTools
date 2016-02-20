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

% Using CKF
filter_opts.use_EKF = 0;
filter_opts.use_SNC = 1;
filter_opts.use_smoother = true;

storage = KalmanFilter(state, meas_store, filter_opts);
plot_cov_err_envelope(storage.cov_store, storage.state_store - true_state*1e3)
title('CKF State Error, with covariance envelope')
state_error = storage.smoothed_state_store - true_state*1e3;
plot_cov_err_envelope(storage.P_smoothed_diag, state_error)
title('CKF Smoothed State Error, with smoothed covariance envelope')

CKF_RMS = sqrt(sum(state_error.*state_error,2)./num_obs);

CKF_pos_err_fig = figure;
CKF_vel_err_fig = figure;
aaxes = {'x', 'y', 'z';'vx', 'vy', 'vz'};
for ii = 1:3
    figure(CKF_pos_err_fig)
    subplot(3,1,ii)
    hold on
    plot(state_error(ii,:))
    plot(storage.P_smoothed_diag(ii,:), 'r--')
    plot(-storage.P_smoothed_diag(ii,:), 'r--')
    title(sprintf('RMS = %.4f m',CKF_RMS(ii)));
    ylabel([aaxes{1,ii} ' err (m)'])
    
    figure(CKF_vel_err_fig)
    subplot(3,1,ii)
    hold on
    plot(state_error(ii+3,:))
    plot(storage.P_smoothed_diag(ii+3,:), 'r--')
    plot(-storage.P_smoothed_diag(ii+3,:), 'r--')
    title(sprintf('RMS = %.4f m/s',CKF_RMS(ii+3)));
    ylabel([aaxes{2,ii} ' err (m)'])
end
figure(CKF_pos_err_fig);
xlabel('Observation')
figure(CKF_vel_err_fig);
xlabel('Observation')
    
batch_store = BatchLeastSquares(state, meas_store, filter_opts);
[~,unique_times,~] = unique(meas_store(:,2));
batch_STM_unique = batch_store.STM_store(:,:,unique_times);
len = length(batch_STM_unique);
batch_state = [];
for idx = 1:len
    batch_state(1:6,idx) = batch_STM_unique(1:6,1:6,idx)*batch_store.x_est(1:6)+ batch_store.ref_state(idx, 1:6)';
end
figure
unique_true_state = true_state(:,unique_times)*1e3;
unique_CKF_STM = storage.STM_accum_store(:,:,unique_times);
plot(arrayfun(@(idx) norm(batch_state(1:6,idx) - unique_true_state(1:6,idx)), 1:len))
figure
plot(arrayfun(@(idx) norm(diag(batch_store.STM_store(:,:,idx))-diag(storage.STM_accum_store(:,:,idx))),1:len))
% pos_errors = arrayfun(@(idx) norm(state_err(1:3,idx)), 1:len);

%% EKF smoothing
filter_opts.use_EKF = 1;
filter_opts.use_SNC = 1;
filter_opts.use_smoother = true;

EKF_storage = KalmanFilter(state, meas_store, filter_opts);
plot_cov_err_envelope(EKF_storage.cov_store, EKF_storage.state_store - true_state*1e3)
plot_cov_err_envelope(EKF_storage.P_smoothed_diag, EKF_storage.x_l_k_store - true_state*1e3)
title('EKF Smoothed State Error, with smoothed covariance envelope')


state_error = EKF_storage.x_l_k_store - true_state*1e3;
EKF_RMS = sqrt(sum(state_error.*state_error,2)./num_obs);

EKF_pos_err_fig = figure;
EKF_vel_err_fig = figure;
for ii = 1:3
    figure(EKF_pos_err_fig)
    subplot(3,1,ii)
    hold on
    plot(state_error(ii,:))
    plot(EKF_storage.P_smoothed_diag(ii,:), 'r--')
    plot(-EKF_storage.P_smoothed_diag(ii,:), 'r--')
    title(sprintf('RMS = %.4f m',EKF_RMS(ii)));
    ylabel([aaxes{1,ii} ' err (m)'])
    
    figure(EKF_vel_err_fig)
    subplot(3,1,ii)
    hold on
    plot(state_error(ii+3,:))
    plot(EKF_storage.P_smoothed_diag(ii+3,:), 'r--')
    plot(-EKF_storage.P_smoothed_diag(ii+3,:), 'r--')
    title(sprintf('RMS = %.4f m/s',EKF_RMS(ii+3)));
    ylabel([aaxes{2,ii} ' err (m)'])
end
figure(EKF_pos_err_fig);
xlabel('Observation')
figure(EKF_vel_err_fig);
xlabel('Observation')
