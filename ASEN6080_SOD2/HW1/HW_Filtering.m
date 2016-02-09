%% Filter Data

%% Initialize
% clearvars -except function_list 
global function_list;
function_list = {};
addpath('C:\Users\John\Documents\Astro\ASEN5070_SOD\tools\')
addpath('C:\Users\John\Documents\Astro\ASEN5050\tools\')

% Load the params file
filter_params;
%     storage = KalmanFilter(state, meas_store, filter_opts);

%% SNC in inertial
SNC_sigma = 1e-1;
cnt = 0;
sig_store = [];
rms_store = [];
pos_RMS = [];
vel_RMS = [];
skip_obs = 100; % observations to skip when computing error RMS.

run_cases = 0;
run_SNC_in_RIC = 0;
if run_cases
while SNC_sigma >= 1e-15
    cnt = cnt+1;
    
    SNC_sigma = SNC_sigma*1e-1;
    sig_store(cnt) = SNC_sigma;
    filter_opts.SNC_Q = eye(3)*SNC_sigma*SNC_sigma;
    storage = KalmanFilter(state, meas_store, filter_opts);
    rms_store(cnt) = storage.RMS;
    state_err = storage.state_store- true_state*1e3;
    len = max(size(state_err));
    pos_errors = arrayfun(@(idx) norm(state_err(1:3,idx)), 1:len);
    vel_errors = arrayfun(@(idx) norm(state_err(4:6,idx)), 1:len);
    pos_errors = pos_errors(skip_obs:end);
    vel_errors = vel_errors(skip_obs:end);
    pos_RMS(cnt) = sqrt(sum(pos_errors.*pos_errors)/len);
    vel_RMS(cnt) = sqrt(sum(vel_errors.*vel_errors)/len);
    fprintf('Finished sigma = %.0e\n', SNC_sigma);
end
end
% SNC in RIC
if run_SNC_in_RIC
    filter_opts.SNC_Q = eye(3)*1e-8;
    filter_opts.SNC_use_RIC = 1;
    storage = KalmanFilter(state, meas_store, filter_opts);
end

% storage = KalmanFilter(state, meas_store, filter_opts);
% plot_cov_err_envelope(storage.cov_store, storage.state_store - true_state*1e3)
% return
    
%% DMC
filter_opts.use_SNC = 0;
filter_opts.use_DMC = 1;
cnt = 1;
% sig_cases = [1e-2 1e-5 1e-7 1e-9 1e-10 1e-11 1e-12 1e-13 1e-15];
sig_cases = [1e-7 1e-7 1e-9 1e-10 1e-11 1e-12 1e-13 1e-15];
for DMC_sig = sig_cases
    cnt = cnt+1;
    filter_opts.DMC.sigma = [1 1 1]*DMC_sig;
    filter_opts.DMC.q_u = diag(filter_opts.DMC.sigma.*filter_opts.DMC.sigma);
    tic;
    storage = KalmanFilter(state, meas_store, filter_opts);
    toc;
    rms_store(cnt) = storage.RMS;
    state_err = storage.state_store- true_state*1e3;
    len = max(size(state_err));
    pos_errors = arrayfun(@(idx) norm(state_err(1:3,idx)), 1:len);
    vel_errors = arrayfun(@(idx) norm(state_err(4:6,idx)), 1:len);
    pos_errors = pos_errors(skip_obs:end);
    vel_errors = vel_errors(skip_obs:end);
    pos_RMS(cnt) = sqrt(sum(pos_errors.*pos_errors)/len);
    vel_RMS(cnt) = sqrt(sum(vel_errors.*vel_errors)/len);
    fprintf('Finished sigma = %.0e\n', DMC_sig);
end
