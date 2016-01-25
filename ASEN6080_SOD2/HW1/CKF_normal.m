%% HW 1 Problem 2: Sequential Processor 
%% Initialize
% clearvars -except function_list 
global function_list;
function_list = {};
% close all

stat_od_proj_init
% ObsData = load('ObsData.txt');
ObsData = meas_store;
obs_r_idx = 3;
obs_rr_idx = 4;
obs_t_idx = 2;
obs_site_idx = 1;
% ObsData(:,obs_r_idx) = ObsData(:,obs_r_idx)*1e3;
% ObsData(:,obs_rr_idx) = ObsData(:,obs_rr_idx)*1e3;

% consts.Re = Re;
% consts.area = drag.A;
% consts.rho = compute_density(ri);
consts.theta_dot = theta_dot;
% consts.m = drag.m;
consts.state_len = 17;

% Initial covariance: no initial correlation between states.
% P0 = eye(consts.state_len)*1e6;
% P0(7,7) = 1e20;
% P0(9:11,9:11) = eye(3)*1e-10;
P0 = eye(consts.state_len)*1e6;
P0(7,7) = 1e10;
P0(9:17,9:17) = eye(9)*1e-10;

% A priori state
% state is composed of [r;v;mu;j2;site1;site2;site3]
x0_ap = zeros(consts.state_len,1);

sig_range = 0.01; % m
sig_rangerate = 0.001; %m/s
W = [1/(sig_range*sig_range) 0; 0 1/(sig_rangerate*sig_rangerate)];
R = [(sig_range*sig_range) 0; 0 (sig_rangerate*sig_rangerate)];

%% Sequential Processor
meas_dt = 10; %sec
times = 0:meas_dt:prop_time;
ode_opts = odeset('RelTol', 1e-12, 'AbsTol', 1e-20);

[num_obs, ~] = size(ObsData);
chol_P0 = chol(P0,'lower');
P0_inv = chol_P0'\inv(chol_P0);
info_mat = P0_inv;
norm_mat = P0_inv*x0_ap;

% Obs. deviation
y1 = zeros(num_obs,1);
y2 = zeros(num_obs,1);

% CKF init
x_est = x0_ap;
P = P0;
obs_time_last = ObsData(1,obs_t_idx);

use_EKF = 0;
EKF_switchover = 100;
use_joseph = 1;
if use_joseph
    P_joseph_store = zeros(num_obs,1);
else
    P_trace_store = zeros(num_obs,1);
end
STM_accum = eye(consts.state_len);
pfr_store = zeros(2,num_obs);
% Run CKF
for ii = 1:num_obs
    site_num = 0;
    for jj = 1:3
        if ObsData(ii, obs_site_idx) == site(jj).id
            site_num = jj;
            break
        end
    end
    
    obs_time = ObsData(ii,obs_t_idx);
    obs_site = ObsData(ii,obs_site_idx);    
    
    % Propagate
    % Get the state on the reference trajectory at this obs time.
    % Get STM from last obs to this one.
    if ii == 1 || obs_time - obs_time_last == 0
        STM_obs2obs = eye(consts.state_len);
        X = state';
    else
        times = obs_time_last:meas_dt:obs_time;
        % Make the STM reflect an epoch time == the last msmnt time
        STM_obs2obs = eye(consts.state_len);
        last_state = [state;...
            reshape(STM_obs2obs(1:important_block(1),1:important_block(2)),...
            important_block(1)*important_block(2),1)];
        
        % Make sure we propagate with the reference values!
        % If a member of the state is in the propagator_opts, set it here.
        % This is especially important for the EKF, since the reference
        % state changes with the estimated deviation.
        propagator_opts.mu = state(7);
        propagator_opts.J2.params.J2 = state(8);
        propagator_opts.J2.params.mu = state(7);
        
        [T,X] = ode45(@two_body_state_dot, times, last_state, ode_opts, ...
            propagator_opts);
        STM_obs2obs(1:important_block(1),1:important_block(2)) = ...
            reshape(X(end,consts.state_len+1:end), ...
            important_block(1), important_block(2));
    end
    state_at_obs = X(end,1:consts.state_len)';
        
    % Calculate measurement deviation y
    t_obs = ObsData(ii,obs_t_idx);
    
    r_comp = compute_range_ECFsite(state_at_obs(1:3),...
        site(site_num).r*1e3,theta_dot*t_obs);
    rr_comp = compute_range_rate_ECFsite(state_at_obs(1:6),...
        site(site_num).r*1e3,theta_dot*t_obs, theta_dot);
    
    y1(ii) = (ObsData(ii,obs_r_idx)-r_comp);
    y2(ii) = (ObsData(ii,obs_rr_idx)-rr_comp);
    
    % Time update
    STM_accum = STM_obs2obs*STM_accum;
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
    H_tilda = stat_od_proj_H_tilda(state_at_obs, consts);
    
    % Kalman gain
    K = P_ap*H_tilda'/(H_tilda*P_ap*H_tilda'+R);
    
    % Measurement Update
    y = [y1(ii);y2(ii)];
    x_est = x_ap + K*(y - H_tilda*x_ap);
    I = eye(consts.state_len);
    if use_joseph
        P = (I-K*H_tilda)*P_ap*(I-K*H_tilda)' + K*R*K';
        P_joseph_store(ii) = trace(P(1:6,1:6));
    else
        P = (I-K*H_tilda)*P_ap;
        P_trace_store(ii) = trace(P(1:6,1:6));
    end
    
    % Track the last time
    obs_time_last = obs_time;    
        
    if ii == 1
        x_est_CKF = x_est;
    end
    pfr = y-H_tilda*x_est;
    pfr_store(:,ii) = pfr;
    
    % Set up ref state for next pass
    state = state_at_obs;
    if use_EKF && ii > EKF_switchover
        state = state + x_est;
        x_est = zeros(consts.state_len,1);
    end

end

RMS_accum = 0;
for ii = 1:num_obs
    RMS_accum = RMS_accum + pfr_store(:,ii)'/R*pfr_store(:,ii);
end
RMS = sqrt(RMS_accum/num_obs)

rangerate_RMS = sqrt(sum(pfr_store(2,:).*pfr_store(2,:))/num_obs)
range_RMS = sqrt(sum(pfr_store(1,:).*pfr_store(1,:))/num_obs)


figure
subplot(2,1,1)
plot(1:num_obs, pfr_store(1,:))
hold on
plot(1:num_obs,3*sig_range*ones(1,num_obs),'r--')
plot(1:num_obs,-3*sig_range*ones(1,num_obs),'r--')
title(sprintf('Range RMS = %.4e m',range_RMS))
ylabel('m')
subplot(2,1,2)
plot(1:num_obs, pfr_store(2,:))
hold on
plot(1:num_obs,3*sig_rangerate*ones(1,num_obs),'r--')
plot(1:num_obs,-3*sig_rangerate*ones(1,num_obs),'r--')
title(sprintf('Range-Rate RMS = %.4e m/s',range_RMS))
ylabel('m/s'),xlabel('Observation')

% %% Covariance Matrix Traces
% % The Joseph formulation follows the same basic shape of the regular Kalman
% % P formulation, but is generally slightly larger. The Joseph formulation
% % will better consider measurements because of this.
% figure
% semilogy(ObsData(:,1)/3600,abs(P_trace_store))
% hold on
% semilogy(ObsData(:,1)/3600,abs(P_joseph_store),'r')
% legend('Kalman P trace', 'Joseph P trace')
% title('Traces of Covariance Matrix')
% xlabel('Time (hr)'), ylabel('S/C Position/Velocity Trace')
% 
% if use_joseph
%     CKF_joseph_init_state = x_est_CKF+state(1:consts.state_len);
%     CKF_joseph_final_state = x_est+X_store(end,1:consts.state_len)';
% else
%     CKF_init_state = x_est_CKF+state(1:consts.state_len);
%     CKF_final_state = x_est+X_store(end,1:consts.state_len)';
%     conv_P_trace_store = P_trace_store;
% end
