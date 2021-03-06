%% HW 1 Problem 2: Sequential Processor 
%% Initialize
% clearvars -except function_list 
% global function_list;
% function_list = {};
% addpath('C:\Users\John\Documents\Astro\ASEN5070_SOD\tools\')
% addpath('C:\Users\John\Documents\Astro\ASEN5050\tools\')
% close all
function output = AttKalmanFilter(state_ap, meas_store, fo)
% filter_params; % load filter params
% ObsData = load('ObsData.txt');
ObsData = meas_store;
obs_t_idx = 2;
hasYaw = false;
if length(ObsData(1,:)) > 4
    hasYaw = true;
end

propagator_opts = fo.propagator_opts;
theta_dot = fo.theta_dot;

consts.theta_dot = fo.theta_dot;
consts.state_len = fo.state_length;

% Initial covariance: no initial correlation between states.
if isfield(fo,'P0')
    P0 = fo.P0;
end

% A priori state
state = state_ap;
PV_state = fo.PV_state;
% Start with zero deviation
x0_ap = zeros(consts.state_len,1);

if isfield(fo,'R')
    R = fo.R;
end

%% Sequential Processor
meas_dt = 10; %sec
fo.ode_opts = odeset('RelTol', 1e-12, 'AbsTol', 1e-20);

[num_obs, ~] = size(ObsData);
% num_obs = 1000; % Just do a subset

% Obs. deviation
y1 = zeros(num_obs,1);
y2 = zeros(num_obs,1);
y3 = zeros(num_obs,1);
if hasYaw
    y4 = zeros(num_obs,1);
end
if fo.YPR_rates
    y4 = zeros(num_obs,1);
    y5 = zeros(num_obs,1);
    y6 = zeros(num_obs,1);
end

% CKF init
x_est = x0_ap;
P = P0;
obs_time_last = ObsData(1,obs_t_idx);

% if fo.use_joseph
%     P_joseph_store = zeros(num_obs,1);
% else
%     P_trace_store = zeros(num_obs,1);
% end
num_state_store = consts.state_len;
if fo.use_EKF

end

% Set up storage
x_est_store = zeros(num_state_store,num_obs);
prefit_range_store = zeros(num_obs,1);
state_store = zeros(num_state_store,num_obs);
state_ap_store = zeros(num_state_store,num_obs);
num_variance_store = num_state_store;
EKF_cov_store = zeros(num_variance_store,num_obs);
if fo.MEKF
    EKF_cov_store = zeros(3,num_obs);
end

% Smoother
% Set up all the storage mechanisms for good memory managements
if fo.use_smoother
    % Store
    % Deviation estimates are already stored
    % STM between observations
    STM_smooth_store = zeros(consts.state_len,consts.state_len,num_obs);
    % Probably a smarter way to store these since they are symmetric.
    P_smooth_store = zeros(consts.state_len,consts.state_len,num_obs);
    P_ap_smooth_store = zeros(consts.state_len,consts.state_len,num_obs);
    x_l_k_store = zeros(consts.state_len,num_obs);
    P_smoothed_diag = zeros(consts.state_len,num_obs);
    STM_accum_store = zeros(consts.state_len,consts.state_len,num_obs);
    b_store = zeros(consts.state_len,num_obs);
    b_last = zeros(consts.state_len,1);
    smoother_propagator_opts = propagator_opts;
    smoother_propagator_opts.OD.use = 0;
end

STM_accum = eye(consts.state_len);
pfr_store = zeros(3,num_obs);
if hasYaw
    pfr_store = zeros(4,num_obs);
end
if fo.YPR_rates
    pfr_store = zeros(6,num_obs);
end
if fo.MEKF
    pfr_store = zeros(3,num_obs);
end

mystr = '';
fprintf('Observation ');

% Run CKF
for ii = 1:num_obs
    fprintf(repmat('\b',size(mystr)))
    mystr = sprintf('%d', ii);
    fprintf(mystr);
    
    obs_time = ObsData(ii,1);
    
    % Propagate
    % Get the state on the reference trajectory at this obs time.
    % Get STM from last obs to this one.
    if fo.MEKF
        if ii == 1 || obs_time - obs_time_last == 0
            STM_obs2obs = eye(consts.state_len);
            X = [state; reshape(P, 9,1)]';
        else
            times = [obs_time_last obs_time];
            last_state = [state;...
                reshape(P, 9,1)];
            [~,X] = ode45(@rigid_body_quat_MEKF, times, last_state, ...
                fo.ode_opts, propagator_opts.att_prop_opts);
        end
        ref_state_at_obs = X(end,1:7)';
        
    else
        if ii == 1 || obs_time - obs_time_last == 0
            STM_obs2obs = eye(consts.state_len);
            X = [PV_state; state]';
        else
    %         times = obs_time_last:meas_dt:obs_time;
            times = [obs_time_last obs_time];
            % Make the STM reflect an epoch time == the last msmnt time
            STM_obs2obs = eye(consts.state_len);
            last_state = [PV_state; state;...
                reshape(STM_obs2obs(1:fo.important_block(1),1:fo.important_block(2)),...
                fo.important_block(1)*fo.important_block(2),1)];

            try
                [~,X] = ode45(@combined_state_dot_estimation, times, last_state, ...
                fo.ode_opts, propagator_opts);
            STM_obs2obs(1:fo.important_block(1),1:fo.important_block(2)) = ...
                reshape(X(end,end-(consts.state_len*consts.state_len-1):end), ...
                fo.important_block(1), fo.important_block(2));
            PV_state = X(end,1:6)';
    %         if fo.use_DMC
    %             STM_obs2obs = compute_DMC_STM(fo, dt, STM_obs2obs(1:fo.important_block(1),1:fo.important_block(2)));
    %         end
            catch
                fprintf('Got that error.\n');
                times
                last_state
                ii
                return
            end
        end
        ref_state_at_obs = X(end,7:7 + consts.state_len-1)';
    end
        
    % Calculate measurement deviation y
    if fo.YPR
        y1(ii) = (ObsData(ii,2)-ref_state_at_obs(1));
        y2(ii) = (ObsData(ii,3)-ref_state_at_obs(2));
        y3(ii) = (ObsData(ii,4)-ref_state_at_obs(3));
        if fo.YPR_rates
            y4(ii) = (ObsData(ii,5)-ref_state_at_obs(4));
            y5(ii) = (ObsData(ii,6)-ref_state_at_obs(5));
            y6(ii) = (ObsData(ii,7)-ref_state_at_obs(6));
        end
    elseif fo.MEKF
        if ref_state_at_obs(1) < 0
            ref_state_at_obs(1:4) = ref_state_at_obs(1:4)*-1;
        end
        dq_z1 = subEP(ObsData(ii,2:5),ref_state_at_obs(1:4));
        dq_z2 = subEP(ObsData(ii,2:5),-ref_state_at_obs(1:4));
%         dq_z = subEP(ref_state_at_obs(1:4),ObsData(ii,2:5));
        if dq_z1(1) > 0
            y1(ii) = 2*dq_z1(2);
            y2(ii) = 2*dq_z1(3);
            y3(ii) = 2*dq_z1(4);   
        else
            y1(ii) = 2*dq_z2(2);
            y2(ii) = 2*dq_z2(3);
            y3(ii) = 2*dq_z2(4);   
        end
%         y1(ii) = 2*(ObsData(ii,2)-ref_state_at_obs(2));
%         y2(ii) = 2*(ObsData(ii,3)-ref_state_at_obs(3));
%         y3(ii) = 2*(ObsData(ii,4)-ref_state_at_obs(4));        
    else
        y1(ii) = (ObsData(ii,2)-ref_state_at_obs(4));
        y2(ii) = (ObsData(ii,3)-ref_state_at_obs(5));
        y3(ii) = (ObsData(ii,4)-ref_state_at_obs(6));
        if hasYaw
            y4(ii) = (ObsData(ii,5)-ref_state_at_obs(1));
        end
    end
        
    
    
    % Time update
    STM_accum = STM_obs2obs*STM_accum;
    if fo.MEKF
        x_ap = zeros(3,1);
        P_ap = reshape(X(end,8:end), 3, 3);
    else
        x_ap = STM_obs2obs*x_est;
        P_ap = STM_obs2obs*P*STM_obs2obs';
    end
    
    dt = obs_time - obs_time_last;
    if fo.use_SNC ...
            && dt < fo.SNC_meas_separation_threshold
        
        % Compute Q. Use the state from last obs for any transformations.
        Q = fo.SNC_Q;
        Gamma = fo.SNC_Gamma(dt);

        P_ap = P_ap + Gamma*Q*Gamma'; % Add the process noise
%     elseif fo.use_DMC
%         Q = compute_DMC_Q(fo,dt);
%         P_ap = P_ap + Q; % DMC process noise
    end
    
    % H~
    consts.t = obs_time;
%     H_tilda = fo.H_tilda_handle(ref_state_at_obs, consts);

    if fo.MEKF
        H_tilda = eye(3);
    elseif fo.YPR
        if fo.YPR_rates
            H_tilda = eye(6);
        else
            H_tilda = [eye(3) zeros(3,consts.state_len-3)];
        end
    else
        H_tilda = [zeros(3,consts.state_len-3) eye(3)];
        if hasYaw
            H_tilda = [H_tilda; 1 0 0 0 0 0];
        end
    end
    
%     if fo.use_DMC
%         H_tilda = [H_tilda zeros(2,3)]; %#ok<AGROW>
%     end
    
    % Kalman gain
    K = P_ap*H_tilda'/(H_tilda*P_ap*H_tilda'+R);
    
    % Measurement Update
    y = [y1(ii);y2(ii);y3(ii)];
    if hasYaw && ~fo.MEKF
        y = [y; y4(ii)];
    end
    if fo.YPR && fo.YPR_rates
        y = [y1(ii);y2(ii);y3(ii);y4(ii);y5(ii);y6(ii)];
    end
    x_est = x_ap + K*(y - H_tilda*x_ap);
    I = eye(consts.state_len);
    if fo.MEKF
        I = eye(3);
    end
    if fo.use_joseph
        P = (I-K*H_tilda)*P_ap*(I-K*H_tilda)' + K*R*K';
%         P_joseph_store(ii) = trace(P(1:6,1:6));
    else
        P = (I-K*H_tilda)*P_ap;
%         P_trace_store(ii) = trace(P(1:6,1:6));
    end
    
    % Track the last time
    obs_time_last = obs_time;    
        
    if ii == 1
        x_est_CKF = x_est;
    end
    pfr = y-H_tilda*x_est;
    pfr_store(:,ii) = pfr;
    
    % Set up ref state for next pass
    if fo.MEKF
        dq = 0.5*[sqrt(4-norm(x_est)^2); x_est];
        state(1:4) = addEP(ref_state_at_obs(1:4),dq);
        state(1:4) = state(1:4)/norm(state(1:4)); %Normalize!
        % Shortest rotation
        if state(1) < 0
            state(1:4) = state(1:4)*-1;
        end
        state(5:7) = ObsData(ii,6:8)';
        
        % Store things
        state_store(:,ii) = state(1:num_state_store);
%         EKF_cov_store(:,ii) = diag(P(1:3,1:3));
        EKF_cov_store(:,ii) = diag(P);
    else
    state = ref_state_at_obs;
    if fo.use_EKF && ii > fo.EKF_switchover
        state = state + x_est;
        x_est_store(:,ii) = state;
        if fo.use_smoother && ii ~= fo.EKF_switchover + 1
            state_ap_store(:,ii) = ref_state_at_obs(1:num_state_store);
        else
            state_ap_store(:,ii) = ref_state_at_obs(1:num_state_store)+x_ap;
        end
        x_est = zeros(consts.state_len,1);
        % Store things
        state_store(:,ii) = state(1:num_state_store);
%         EKF_cov_store(:,ii) = diag(P(1:3,1:3));
        EKF_cov_store(:,ii) = diag(P);
    else%if ~fo.use_EKF && ii > fo.EKF_switchover
        % Nothing to calculate, just store.
        x_est_store(:,ii) = x_est;
        if fo.use_smoother && fo.use_EKF
            state_ap_store(:,ii) = ref_state_at_obs(1:num_state_store)+x_ap;
        end
        prefit_range_store(ii) = y1(ii);
        state_store(:,ii) = ...
            state(1:num_state_store) + x_est(1:num_state_store);
        EKF_cov_store(:,ii) = ...
            diag(P(1:num_variance_store,1:num_variance_store));
    end
    end
    
    if fo.use_smoother
        % Store
        STM_smooth_store(:,:,ii) = STM_obs2obs;
        P_smooth_store(:,:,ii) = P;
        P_ap_smooth_store(:,:,ii) = P_ap;
        STM_accum_store(:,:,ii) = STM_accum;
    end
end
fprintf('\n');

RMS_accum = 0;
for ii = 1:num_obs
    RMS_accum = RMS_accum + pfr_store(:,ii)'/R*pfr_store(:,ii);
end
output.RMS = sqrt(RMS_accum/num_obs)

output.rangerate_RMS = sqrt(sum(pfr_store(2,:).*pfr_store(2,:))/num_obs)
output.range_RMS = sqrt(sum(pfr_store(1,:).*pfr_store(1,:))/num_obs)
output.pfr_store = pfr_store;
output.prefit_range_store = prefit_range_store;
output.cov_store = EKF_cov_store;
output.state_store = state_store;
output.x_est_store = x_est_store;
output.state_ap_store = state_ap_store;
output.final_P = P;

%% Smoothed state
if fo.use_smoother
    % deviation at time k with information up to l (the end)
    x_l_k = x_est_store(:,end);
    if fo.use_EKF
        x_l_k = state_store(:,end);
    end
    x_l_k_store(:,end) = x_l_k;
    P_l_k = P_smooth_store(:,:,end);
    P_smoothed_diag(:,end) = diag(P_l_k);
    if fo.use_EKF
        bk = zeros(consts.state_len,1);
    end
    for ii = 1:(num_obs-1) 
        Sk = P_smooth_store(:,:,end-ii) * STM_smooth_store(:,:,end-ii+1)' ...
            /P_ap_smooth_store(:,:,end-ii+1);
        if fo.use_EKF
            x_l_k = state_store(:,end-ii) ...
                    + Sk*(x_l_k - state_ap_store(:,end-ii+1));
        else
            x_l_k = x_est_store(:,end-ii) ...
                + Sk*(x_l_k - STM_smooth_store(:,:,end-ii+1)*x_est_store(:,end-ii));
            
        end
        x_l_k_store(:,end-ii) = x_l_k;
        P_l_k = P_smooth_store(:,:,end-ii) + Sk*(P_l_k - P_ap_smooth_store(:,:,end-ii+1))*Sk';
        P_smoothed_diag(:,end-ii) = diag(P_l_k);
    end
    if fo.use_EKF
        output.b = b_store;
    end
    smoothed_state_store = state_store - x_est_store + x_l_k_store;
    output.x_l_k_store = x_l_k_store;
    output.smoothed_state_store = smoothed_state_store;
    output.P_smoothed_diag = P_smoothed_diag;
    output.STM_smooth_store = STM_smooth_store;
    output.STM_accum_store = STM_accum_store;
end

%%
% sig_range = 1e-6;
% sig_rangerate = 1e-6;
% figure
% subplot(2,1,1)
% plot(1:num_obs, pfr_store(1,:),'.','LineWidth',1)
% hold on
% plot(1:num_obs,3*sig_range*ones(1,num_obs),'r--')
% plot(1:num_obs,-3*sig_range*ones(1,num_obs),'r--')
% title(sprintf('Range RMS = %.4e m',output.range_RMS))
% ylabel('m')
% subplot(2,1,2)
% plot(1:num_obs, pfr_store(2,:),'.','LineWidth',1)
% hold on
% plot(1:num_obs,3*sig_rangerate*ones(1,num_obs),'r--')
% plot(1:num_obs,-3*sig_rangerate*ones(1,num_obs),'r--')
% title(sprintf('Range-Rate RMS = %.4e m/s',output.rangerate_RMS))
% ylabel('m/s'),xlabel('Observation')

% %% Covariance Matrix Traces
% % The Joseph formulation follows the same basic shape of the regular Kalman
% % P formulation, but is generally slightly larger. The Joseph formulation
% % will better consider measurements because of this.
% figure
% semilogy(ObsData(:,1)/3600,abs(P_trace_store))
% hold on
% semilogy(ObsData(:,2)/3600,abs(P_joseph_store),'r')
% legend('Kalman P trace', 'Joseph P trace')
% title('Traces of Covariance Matrix')
% xlabel('Time (hr)'), ylabel('S/C Position/Velocity Trace')
% 
% if fo.use_joseph
%     CKF_joseph_init_state = x_est_CKF+state(1:consts.state_len);
%     CKF_joseph_final_state = x_est+X_store(end,1:consts.state_len)';
% else
%     CKF_init_state = x_est_CKF+state(1:consts.state_len);
%     CKF_final_state = x_est+X_store(end,1:consts.state_len)';
%     conv_P_trace_store = P_trace_store;
% end

% Sound when done
% [y_sound,Fs] = audioread('C:\Users\John\Documents\StarCraft_Sound_Pack\Protoss\Units\Artanis\patpss00.wav');
% sound(y_sound,Fs);
end

%% SNC process noise
function Q = compute_SNC_Q(fo, X)
    if fo.SNC_use_RIC
        % Find the rotation matrix RIC->Inertial
        r_inrtl = X(1:3)/norm(X(1:3));
        h = cross(X(1:3),X(4:6));
        cross_track_inrtl = h/norm(h);
        in_track = cross(cross_track_inrtl, r_inrtl);
        
        R_eci_ric = [r_inrtl in_track cross_track_inrtl];
        Q = R_eci_ric'*fo.SNC_Q*R_eci_ric;
    else
        Q = fo.SNC_Q;
    end

end

%% DMC process noise
function Q = compute_DMC_Q(fo,dt)
    % Set up some vars for faster computation
    beta = diag(fo.DMC.B);
    beta_2 = beta.*beta;
    beta_3 = beta.*beta_2;
    beta_4 = beta.*beta_3;
    beta_5 = beta.*beta_4;
    sig_sq = diag(fo.DMC.q_u);
    dt_2 = dt*dt;
    dt_3 = dt*dt_2;
    exp_term = exp(-beta*dt);
    exp_sq_term = exp(-2*beta*dt);
    
    % Fill out all the sub-matrices
    Q_rr = diag(sig_sq.*...
        (dt_3./(3.*beta_2) - dt_2./beta_3 ...
        + dt./beta_4 - 2*exp_term*dt./beta_4 ...
        + (1-exp_sq_term)./(2.*beta_5)));
    
    % Q_rv = Q_vr
    Q_rv = diag(sig_sq.*...
        (dt_2./(2.*beta_2) - dt./beta_3 + exp_term*dt./beta_3 ...
        +(1-exp_term)./beta_4 - (1-exp_sq_term)./(2.*beta_4)));
    
    %Q_rw = Q_wr
    Q_rw = diag(sig_sq.*...
        ((1-exp_sq_term)./(2.*beta_3) - exp_term*dt./beta_2));
    
    Q_vv = diag(sig_sq.*...
        (dt./beta_2 ...
        - 2*(1-exp_term)./beta_3 + (1-exp_sq_term)./(2*beta_3)));
    
    %Q_vw = Q_wv
    Q_vw = diag(sig_sq.*...
        ((1+exp_sq_term)./(2*beta_2) - exp_term./beta_2));
    
    Q_ww = diag(sig_sq.*(1-exp_sq_term)./(2*beta));
    
    Q = [Q_rr Q_rv Q_rw;
         Q_rv Q_vv Q_vw;
         Q_rw Q_vw Q_ww];
end

%% DMC STM calculation.
function STM = compute_DMC_STM(fo, dt, STM_non_aug)
    % Set up some vars for faster computation
    beta = diag(fo.DMC.B);
    beta_2 = beta.*beta;
    exp_term = exp(-beta*dt);
    M = diag(exp_term);
    PhiWM = [diag((1-exp_term)./beta); ...
        diag((exp_term-1)./beta_2 + dt./beta)];
    STM = [STM_non_aug PhiWM; zeros(3,6) M];
        
end