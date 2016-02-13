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
filter_opts.use_EKF = 1;

storage = KalmanFilter(state, meas_store, filter_opts);
plot_cov_err_envelope(storage.cov_store, storage.state_store - true_state*1e3)
% return
    
