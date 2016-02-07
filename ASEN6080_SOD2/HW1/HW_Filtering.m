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

run_cases = 0;
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
    pos_RMS(cnt) = sqrt(sum(pos_errors.*pos_errors)/len);
    vel_RMS(cnt) = sqrt(sum(vel_errors.*vel_errors)/len);
    fprintf('Finished sigma = %.0e\n', SNC_sigma);
end
end

filter_opts.SNC_Q = eye(3)*1e-8;
filter_opts.SNC_use_RIC = 1;
storage = KalmanFilter(state, meas_store, filter_opts);

