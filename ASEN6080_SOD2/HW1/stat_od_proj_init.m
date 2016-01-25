fcnPrintQueue(mfilename('fullpath'))
%% Initial data for stat OD project
mu = 3.986004415e14; %m3/s2
J2 = 1.082626925638815e-3;
J3 =  -2.5327e-6;
Re = 6378136.3; %m

theta_dot = 7.2921158553e-5; %rad/s
% 
% site(1).name = 'Pacific Ocean Ship Sensor';
% site(1).id = 101;
% site(1).r = [-5127510.0 -3794160.0 0.0]'; % m
% 
% site(2).name = 'Pirinclik, Turkey';
% site(2).id = 337;
% site(2).r = [3860910.0 3238490.0 3898094.0]'; % m
% 
% site(3).name = 'Thule, Greenland';
% site(3).id = 394;
% site(3).r = [549505.0 -1380872.0 6182197.0]'; % m

% ri = [757700.0 5222607.0 4851500.0]';
% vi = [2213.21 4678.34 -5371.30]';

% x0_ap(1:3) = x0_ap(1:3) - 1000; %m
% x0_ap(4:6) = x0_ap(4:6) + 15; %m/s
drag.Cd = 0.0;
drag.A = 3.0; % m
drag.m = 970; %kg

% state = [ri; vi; mu; J2; drag.Cd; site(1).r; site(2).r; site(3).r];
state = [state_i*1e3; mu; J2; site(1).r*1e3; site(2).r*1e3; site(3).r*1e3];

% Set up propagator options
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
propagator_opts.OD.state_len = 17;
propagator_opts.OD.A_mat_handle = @stat_od_proj_A;
propagator_opts.OD.A_params.Re = Re;
propagator_opts.OD.A_params.area = drag.A;
propagator_opts.OD.A_params.rho = propagator_opts.drag.model_params.rho0;
propagator_opts.OD.A_params.theta_dot = theta_dot;
propagator_opts.OD.A_params.m = drag.m;
propagator_opts.OD.A_params.H = propagator_opts.drag.model_params.H;
important_block = [8 8]; %rows, cols
propagator_opts.OD.A_params.important_block = important_block;
propagator_opts.OD.A_params.state_len = propagator_opts.OD.state_len;
STM_i = eye(propagator_opts.OD.state_len);
% state = [state; reshape(STM_i(1:important_block(1),1:important_block(2)),...
%     important_block(1)*important_block(2),1)];