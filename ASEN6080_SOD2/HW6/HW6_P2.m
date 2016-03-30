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

filter_opts.use_EKF = 0;
filter_opts.use_SNC = 0;
filter_opts.use_DMC = 0;
filter_opts.use_joseph = 0;

tic
CKF_storage = ConsiderKalmanFilter(state, meas_store, filter_opts);
toc
plot_cov_err_envelope(CKF_storage.cov_store(1:3,:), CKF_storage.state_store - true_state*1e3)
plot_cov_err_envelope(CKF_storage.consider_cov_store(1:3,:), CKF_storage.state_store - true_state*1e3)
residual_plot(CKF_storage.pfr_store, sigs,'CKF')

mapped_x_est0 = CKF_storage.STM\CKF_storage.x_est_store(:,end);

state_0 = state + mapped_x_est0;

propagator_opts.OD.use = 0;
[C, IA, IC] = unique(meas_store(:,2));
[~,X] = ode45(@two_body_state_dot, C, state_0, ...
            odeset('RelTol', 1e-12, 'AbsTol', 1e-20), propagator_opts);
full_Pc0 = ...\
    CKF_storage.Psi(:,:,end)\CKF_storage.final_full_Pc/(CKF_storage.Psi(:,:,end)');
prop_Pc_store = zeros(7,num_obs);
for ii = 1:num_obs
    Pck =  CKF_storage.Psi(:,:,ii)*full_Pc0*(CKF_storage.Psi(:,:,ii)');
    prop_Pc_store(:,ii) = diag(Pck);
end
plot_cov_err_envelope(prop_Pc_store(1:3,:), CKF_storage.state_store - [X(IC,1) X(IC,2) X(IC,3) X(IC,4) X(IC,5) X(IC,6)]')
        