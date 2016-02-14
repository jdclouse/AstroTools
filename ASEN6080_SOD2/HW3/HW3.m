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
skip_obs = 100; % observations to skip when computing error RMS.

% Using CKF
filter_opts.use_EKF = 0;
filter_opts.use_SNC = 1;
filter_opts.use_smoother = true;

storage = KalmanFilter(state, meas_store, filter_opts);
plot_cov_err_envelope(storage.cov_store, storage.state_store - true_state*1e3)
plot_cov_err_envelope(storage.P_smoothed_diag, storage.smoothed_state_store - true_state*1e3)
% return
    
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