%% HW5

%% Initialize
% clearvars -except function_list 
global function_list;
function_list = {};
addpath('C:\Users\John\Documents\Astro\ASEN6080_SOD\tools\')

% Load the params file
filter_params;

%% 
num_obs = length(meas_store);
sig_range = 0.01; % m
sig_rangerate = 0.001; %m/s
sigs = [sig_range; sig_rangerate];

% Using UKF
% P = eye(6)*1e2;
% filter_opts.use_SNC = 0;
filter_opts.propagator_opts.J3.use = 0;
tic
storage = SRIF(state, P, meas_store, filter_opts);
toc

