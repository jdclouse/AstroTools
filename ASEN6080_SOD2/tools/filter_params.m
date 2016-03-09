fcnPrintQueue(mfilename('fullpath'))
%% Initial data for stat OD project
% mu = 3.986004415e14; %m3/s2
mu = 3.986004418e14; %m3/s2
% J2 = 1.082626925638815e-3;
J2 = 0.0010826267;
J3 =  -2.5327e-6;
Re = 6378136.3; %m

filter_opts.theta_dot = 7.2921158553e-5; %rad/s % TODO clean this!
filter_opts.num_sites = num_sites;
filter_opts.site = site;

drag.Cd = 0.0;
drag.A = 3.0; % m
drag.m = 970; %kg

state = state_i*1e3;
state(1:3) = state(1:3) + 50;
state(4:6) = state(4:6) - 0.01;
P = eye(6)*1e6;

% Set up propagator options for the filter
propagator_opts.param_in_state.mu_idx = 0;
propagator_opts.param_in_state.J2_idx = 0;
propagator_opts.param_in_state.J3_idx = 0;
propagator_opts.param_in_state.Cd_idx = 0;
propagator_opts.mu = mu; 

propagator_opts.drag = drag;
propagator_opts.drag.use = 0;
propagator_opts.drag.model_params.rho0 = 3.614e-13; %kg/m3
propagator_opts.drag.model_params.r0 = 700000+6378136.3;
propagator_opts.drag.model_params.H = 88667;
propagator_opts.drag.model_params.theta_dot = theta_dot;

propagator_opts.J2.use = 1;
propagator_opts.J2.params.J2 = J2;
propagator_opts.J2.params.mu = mu; 
propagator_opts.J2.params.Re = Re;
propagator_opts.J3.use = 0;
propagator_opts.J3.params.J3 = J3;
propagator_opts.J3.params.mu = mu; 
propagator_opts.J3.params.Re = Re;

propagator_opts.OD.use = 1;
propagator_opts.OD.state_len = 6;
propagator_opts.OD.A_mat_handle = @A_state_rv;
propagator_opts.OD.A_params.mu = propagator_opts.mu;
propagator_opts.OD.A_params.J2 = propagator_opts.J2;
propagator_opts.OD.A_params.J3 = propagator_opts.J3;
propagator_opts.OD.A_params.Re = Re;
propagator_opts.OD.A_params.area = drag.A;
propagator_opts.OD.A_params.rho = propagator_opts.drag.model_params.rho0;
propagator_opts.OD.A_params.theta_dot = theta_dot;
propagator_opts.OD.A_params.m = drag.m;
propagator_opts.OD.A_params.H = propagator_opts.drag.model_params.H;
filter_opts.important_block = [6 6]; %rows, cols
propagator_opts.OD.A_params.important_block = filter_opts.important_block;
propagator_opts.OD.A_params.state_len = propagator_opts.OD.state_len;
STM_i = eye(propagator_opts.OD.state_len);
% state = [state; reshape(STM_i(1:important_block(1),1:important_block(2)),...
%     important_block(1)*important_block(2),1)];

% Filter Options
filter_opts.propagator_opts = propagator_opts;
filter_opts.use_EKF = 1;
filter_opts.use_SNC = 1;
filter_opts.EKF_switchover = 200;
filter_opts.use_joseph = 1;
filter_opts.use_potter = 0;
filter_opts.H_tilda_handle = @H_tilda_state_rv;
filter_opts.SNC_Q = eye(3)*1e-6;
filter_opts.SNC_Gamma = @(dt) [dt*dt/2 0 0;...
            0 dt*dt/2 0;...
            0 0 dt*dt/2;...
            dt 0 0;...
            0 dt 0;...
            0 0 dt];
filter_opts.SNC_meas_separation_threshold = 3600;
filter_opts.SNC_use_RIC = 0;
% Eventually add storage options for analysis

% filter_opts.use_DMC = 0;
% filter_opts.DMC.tau = [1 1 1]'*3/orb_period; % Actually this is 1/tau
% filter_opts.DMC.B = diag(filter_opts.DMC.tau);
% filter_opts.DMC.sigma = [1 1 1]*1e-7;
% filter_opts.DMC.q_u = diag(filter_opts.DMC.sigma.*filter_opts.DMC.sigma);
% filter_opts.DMC.w_P0 = eye(3)*1e-5; % a priori covariance

% Smoothing
filter_opts.use_smoother = false;

% UKF
filter_opts.propagator_opts.UKF.L = 6;
filter_opts.UKF.L = 6; 
filter_opts.UKF.alpha = 1; 
filter_opts.UKF.beta = 2; 