%% HW 10: Sequential Processor for the Term Project
%% Initialize
% clearvars -except function_list pub_opt
global function_list;
function_list = {};
close all

stat_od_proj_init
ObsData = load('ObsData.txt');

% Calculate PHI(18340,0)

consts.Re = Re;
consts.area = drag.A;
consts.rho = compute_density(ri);
consts.theta_dot = theta_dot;
consts.m = drag.m;
consts.state_len = 18;

P0 = eye(consts.state_len)*1e6;
P0(7,7) = 1e20;
P0(10:12,10:12) = eye(3)*1e-10;

x0_ap = zeros(consts.state_len,1);

sig_range = 0.01; % m
sig_rangerate = 0.001; %m/s
W = [1/(sig_range*sig_range) 0; 0 1/(sig_rangerate*sig_rangerate)];
R = [(sig_range*sig_range) 0; 0 (sig_rangerate*sig_rangerate)];


dt = 0.1;
times = 0:dt:18340;
ode_opts = odeset('RelTol', 1e-12, 'AbsTol', 1e-20);

% for iter = 1:3
[T,X] = ode45(@two_body_state_dot, times, state, ode_opts, propagator_opts);

% Store off every 20 seconds of data
X_store = X(mod(times,20) == 0,:);
T_store = T(mod(times,20) == 0);

% % State at t = 18340 seconds
% times2 = 
% [T,X] = ode45(@two_body_state_dot, times, X_store(end,:)', ode_opts, propagator_opts);

% Accumulate the information and normal matrices
[num_obs, ~] = size(ObsData);
% info_mat = zeros(consts.state_len);
chol_P0 = chol(P0,'lower');
P0_inv = chol_P0'\inv(chol_P0);
info_mat = P0_inv;
norm_mat = P0_inv*x0_ap;
H_tilda_given = load('BatchHtilda.mat');
cntr =1 ;
% Obs. deviation
y1 = zeros(num_obs,1);
y2 = zeros(num_obs,1);
for ii = 1:num_obs
    site_num = 0;
    for jj = 1:3
        if ObsData(ii, 2) == site(jj).id
            site_num = jj;
            break
        end
    end
    t_obs = ObsData(ii,1);
    ostate = X(T(:,1)==t_obs,1:6);
    
    r_comp = compute_range_ECFsite(ostate(1:3),...
        site(site_num).r,theta_dot*t_obs);
    rr_comp = compute_range_rate_ECFsite(ostate(1:6),...
        site(site_num).r,theta_dot*t_obs, theta_dot);
    
    y1(ii) = (ObsData(ii,3)-r_comp);
    y2(ii) = (ObsData(ii,4)-rr_comp);

end

x_est = x0_ap;
P = P0;
obs_time_last = ObsData(ii,1);
for ii = 1:num_obs
    obs_time = ObsData(ii,1);
    obs_site = ObsData(ii,2);    
    
    % STM from last obs to this one.
    % Not very efficient, since I'm running the integrator again.
    if ii == 1
        STM_obs2obs = eye(consts.state_len);
    else
        times_temp = obs_time_last:dt:obs_time;
        last_state = X_store(T_store == ObsData(ii-1,1),:)';
        STM_obs2obs = eye(consts.state_len);
        % Make the STM reflect an epoch time == the last msmnt time
        last_state(consts.state_len+1:end) = ...
            reshape(STM_obs2obs(1:important_block(1),1:important_block(2)),...
            important_block(1)*important_block(2),1);
        
        [T_temp,X_temp] = ...
            ode45(@two_body_state_dot, times_temp, last_state, ...
            ode_opts, propagator_opts);
        STM_obs2obs(1:important_block(1),1:important_block(2)) = ...
            reshape(X_temp(end,consts.state_len+1:end), ...
            important_block(1), important_block(2));
    end
    obs_time_last = obs_time;    
    
    % Time update
    x_ap = STM_obs2obs*x_est;
    P_ap = STM_obs2obs*P*STM_obs2obs';
    
    % H~
    consts.t = obs_time;
    for xx = 1:3
        if site(xx).id == obs_site
            consts.site = xx;
            break
        end
    end
    state_at_obs = X_store(T_store == obs_time,1:consts.state_len);
    H_tilda = stat_od_proj_H_tilda(state_at_obs, consts);
    
    % Kalman gain
    K = P_ap*H_tilda'/(H_tilda*P_ap*H_tilda'+R);
    
    % Measurement Update
    y = [y1(ii);y2(ii)];
    x_est = x_ap + K*(y - H_tilda*x_ap);
    P = (eye(consts.state_len)-K*H_tilda)*P_ap;
end
STM_18340_0 = eye(consts.state_len);
STM_18340_0(1:important_block(1),1:important_block(2)) = ...
    reshape(X_store(end,consts.state_len+1:end), ...
    important_block(1), important_block(2));
est_x0 = STM_18340_0\x_est;
% x_est = cholesky_linear_solver(info_mat,norm_mat);
% % x0_ap = x0_ap-x_est;
% % state(1:18) = state(1:18) + x_est;
% % end
% 
% info_mat_given = load('BatchInfoMat.mat');
% info_diff = abs((info_mat-info_mat_given.InfoMat)./info_mat_given.InfoMat);
% hist(reshape(log10(info_diff),18*18,1))
% title('Information matrix diff histogram.')
% xlabel('Exponent')
% ylabel('Num Elements')
% 
% figure
% norm_mat_given = load('BatchNormMat.mat');
% norm_diff = abs((norm_mat-norm_mat_given.NormMat)./norm_mat_given.NormMat);
% hist(log10(norm_diff))
% title('Normal matrix diff histogram.')
% xlabel('Exponent')
% ylabel('Num Elements')
% 
x_hat_given = load('BatchXhat');
x_diff = abs((est_x0 - x_hat_given.xhat_pass1)./x_hat_given.xhat_pass1);
est_x0
fprintf('Estimated State Deviation Error\n')
for ii = 1:consts.state_len
    fprintf('%.0e\n',x_diff(ii))
end
