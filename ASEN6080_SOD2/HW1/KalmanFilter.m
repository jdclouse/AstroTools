%% HW 1 Problem 2: Sequential Processor 
%% Initialize
% clearvars -except function_list 
% global function_list;
% function_list = {};
% addpath('C:\Users\John\Documents\Astro\ASEN5070_SOD\tools\')
% addpath('C:\Users\John\Documents\Astro\ASEN5050\tools\')
% close all
function output = KalmanFilter(state_ap, meas_store, fo)
% filter_params; % load filter params
% ObsData = load('ObsData.txt');
ObsData = meas_store;
obs_r_idx = 3;
obs_rr_idx = 4;
obs_t_idx = 2;
obs_site_idx = 1;

propagator_opts = fo.propagator_opts;
site = fo.site;
theta_dot = fo.theta_dot;

consts.theta_dot = fo.theta_dot;
consts.state_len = propagator_opts.OD.state_len;

% Initial covariance: no initial correlation between states.
P0 = eye(consts.state_len)*1e8;
W0 = sqrt(P0);

% A priori state
state = state_ap;
% Start with zero deviation
x0_ap = zeros(consts.state_len,1);

sig_range = 0.01; % m
sig_rangerate = 0.001; %m/s
% W = [1/(sig_range*sig_range) 0; 0 1/(sig_rangerate*sig_rangerate)];
R = [(sig_range*sig_range) 0; 0 (sig_rangerate*sig_rangerate)];

% DMC setup
if fo.use_DMC
    B = diag(fo.DMC.tau);
    C = [zeros(6,3);eye(3)];
    D = [zeros(3);eye(3)];
%     q_u;
    w = [0 0 0]';
    state = [state; w];
    x0_ap = zeros(9,1);
    P0 = [P0 zeros(6,3); zeros(3,6) fo.DMC.w_P0]; %tried q but that doesn't seem like a good idea now.
    consts.state_len = 9;
    fo.important_block = fo.important_block + 3;
    
    propagator_opts.OD.aug_len = 3;
    propagator_opts.OD.DMC.use = 1;
    propagator_opts.OD.DMC.B = fo.DMC.B;
    propagator_opts.OD.DMC.C = C;
    propagator_opts.OD.DMC.D = D;
end

%% Sequential Processor
meas_dt = 10; %sec
fo.ode_opts = odeset('RelTol', 1e-12, 'AbsTol', 1e-20);

[num_obs, ~] = size(ObsData);
% num_obs = 1000; % Just do a subset

% Obs. deviation
y1 = zeros(num_obs,1);
y2 = zeros(num_obs,1);

% CKF init
x_est = x0_ap;
P = P0;
W = W0;
obs_time_last = ObsData(1,obs_t_idx);

if fo.use_joseph
    P_joseph_store = zeros(num_obs,1);
else
    P_trace_store = zeros(num_obs,1);
end
if fo.use_EKF
    EKF_x_est_store = zeros(num_obs,1);
    EKF_prefit_range_store = zeros(num_obs,1);
    EKF_postfit_range_store = zeros(num_obs,1);
    EKF_state_store = zeros(6,num_obs);
else
    CKF_x_est_store = zeros(num_obs,1);
    CKF_prefit_range_store = zeros(num_obs,1);
    CKF_postfit_range_store = zeros(num_obs,1);
    CKF_state_store = zeros(6,num_obs);
end
STM_accum = eye(consts.state_len);
pfr_store = zeros(2,num_obs);

mystr = '';
fprintf('Observation ');

% Run CKF
for ii = 1:num_obs
    fprintf(repmat('\b',size(mystr)))
    mystr = sprintf('%d', ii);
    fprintf(mystr);
    site_num = 0;
    for jj = 1:fo.num_sites
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
            reshape(STM_obs2obs(1:fo.important_block(1),1:fo.important_block(2)),...
            fo.important_block(1)*fo.important_block(2),1)];
        
        % Make sure we propagate with the reference values!
        % If a member of the state is in the propagator_opts, set it here.
        % This is especially important for the EKF, since the reference
        % state changes with the estimated deviation.
        % Mu
        if propagator_opts.param_in_state.mu_idx > 0
            propagator_opts.mu = ...
                state(propatagor_opts.param_in_state.mu_idx);
        end
        % J2
        if propagator_opts.param_in_state.J2_idx > 0
            propagator_opts.J2.params.J2 = ...
                state(propatagor_opts.param_in_state.J2_idx);
            if propatagor_opts.param_in_state.mu_idx > 0
                propagator_opts.J2.params.mu = ...
                    state(propatagor_opts.param_in_state.mu_idx);
            end
        end
        % J3
        if propagator_opts.param_in_state.J3_idx > 0
            propagator_opts.J3.params.J3 = ...
                state(propatagor_opts.param_in_state.J3_idx);
            if propatagor_opts.param_in_state.mu_idx > 0
                propagator_opts.J3.params.mu = ...
                    state(propatagor_opts.param_in_state.mu_idx);
            end
        end
        % Cd
        if propagator_opts.param_in_state.Cd_idx > 0
            propagator_opts.drag.Cd = ...
                state(propatagor_opts.param_in_state.Cd_idx);
        end
        
        try
            [~,X] = ode45(@two_body_state_dot, times, last_state, ...
            fo.ode_opts, propagator_opts);
        STM_obs2obs(1:fo.important_block(1),1:fo.important_block(2)) = ...
            reshape(X(end,consts.state_len+1:end), ...
            fo.important_block(1), fo.important_block(2));
        catch
            fprintf('Got that error.\n');
            times
            last_state
            ii
            return
        end
    end
    ref_state_at_obs = X(end,1:consts.state_len)';
        
    % Calculate measurement deviation y
    t_obs = ObsData(ii,obs_t_idx);
    
    r_comp = compute_range_ECFsite(ref_state_at_obs(1:3),...
        site(site_num).r*1e3,theta_dot*t_obs);
    rr_comp = compute_range_rate_ECFsite(ref_state_at_obs(1:6),...
        site(site_num).r*1e3,theta_dot*t_obs, theta_dot);
    
    y1(ii) = (ObsData(ii,obs_r_idx)-r_comp);
    y2(ii) = (ObsData(ii,obs_rr_idx)-rr_comp);
    
    % Time update
    STM_accum = STM_obs2obs*STM_accum;
    x_ap = STM_obs2obs*x_est;
    P_ap = STM_obs2obs*P*STM_obs2obs';
    
    dt = obs_time - obs_time_last;
    if fo.use_SNC ...
            && dt < fo.SNC_meas_separation_threshold
        
        % Compute Q. Use the state from last obs for any transformations.
        Q = compute_SNC_Q(fo, X(1,1:consts.state_len)');
        Gamma = fo.SNC_Gamma(dt);

        P_ap = P_ap + Gamma*Q*Gamma'; % Add the process noise
    elseif fo.use_DMC
        Q = compute_DMC_Q(fo,dt);
        P_ap = P_ap + Q; % DMC process noise
    end
    
    % H~
    consts.t = obs_time;
    for xx = 1:fo.num_sites
        if site(xx).id == obs_site
            consts.site = xx;
            consts.site_r = site(xx).r*1e3;
            break
        end
    end
    H_tilda = fo.H_tilda_handle(ref_state_at_obs, consts);
    if fo.use_DMC
        H_tilda = [H_tilda zeros(2,3)]; %#ok<AGROW>
    end
    
    % Kalman gain
    K = P_ap*H_tilda'/(H_tilda*P_ap*H_tilda'+R);
    
    % Measurement Update
    y = [y1(ii);y2(ii)];
    x_est = x_ap + K*(y - H_tilda*x_ap);
    I = eye(consts.state_len);
    if fo.use_joseph
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
    state = ref_state_at_obs;
    if fo.use_EKF && ii > fo.EKF_switchover
        EKF_x_est_store(ii) = norm(x_est(1:3));
        state = state + x_est;
        x_est = zeros(consts.state_len,1);
        EKF_state_store(:,ii) = state(1:6);
        EKF_cov_store(:,ii) = diag(P(1:3,1:3));
    elseif ~fo.use_EKF && ii > fo.EKF_switchover
        state;
        CKF_x_est_store(ii) = norm(x_est(1:3));
        CKF_prefit_range_store(ii) = y1(ii);
        CKF_state_store(:,ii) = state(1:6) + x_est(1:6);
        CKF_cov_store(:,ii) = diag(P(1:3,1:3));
    else
        EKF_state_store(:,ii) = state(1:6) + x_est(1:6);
        CKF_state_store(:,ii) = state(1:6) + x_est(1:6);
        EKF_cov_store(:,ii) = diag(P(1:3,1:3));
        CKF_cov_store(:,ii) = diag(P(1:3,1:3));
    end
    if ii > num_obs/2
        state;
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
output.prefit_range_store = EKF_prefit_range_store;
output.cov_store = EKF_cov_store;
output.state_store = EKF_state_store;

% if fo.use_EKF
%     EKF_postfit_range_store = pfr_store;
%     EKF_prefit_range_store = y1;
%     for diff_idx = 1:num_obs
%         cov(diff_idx) = 3*norm(EKF_cov_store(:,diff_idx));
%         diff(diff_idx) = norm(EKF_state_store(:,diff_idx)) - norm(true_state(:,diff_idx)*1e3);
%     end
%     figure; plot(diff);
%     hold on
%     plot(cov,'r')
%     plot(-cov,'r')
%     title('EKF State error, with covariance envelope')
%     xlabel('m')
%     ylabel('Observation')
% else
%     CKF_postfit_range_store = pfr_store;
%     CKF_prefit_range_store = y1;
%     for diff_idx = 1:num_obs
%         cov(diff_idx) = 3*norm(CKF_cov_store(:,diff_idx));
%         diff(diff_idx) = norm(CKF_state_store(:,diff_idx)) - norm(true_state(:,diff_idx)*1e3);
%     end
%     figure; plot(diff);
%     hold on
%     plot(cov,'r')
%     plot(-cov,'r')
%     title('CKF State error, with covariance envelope')
%     xlabel('m')
%     ylabel('Observation')
% end

% figure
%     for diff_idx = 1:num_obs
%         diff(diff_idx) = norm(CKF_state_store(:,diff_idx)) - norm(EKF_state_store(:,diff_idx));
%     end
% plot(diff)
% title('Error between CKF and EKF')
% xlabel('Observation')
% ylabel('m')
%%

figure
subplot(2,1,1)
plot(1:num_obs, pfr_store(1,:),'.','LineWidth',1)
hold on
plot(1:num_obs,3*sig_range*ones(1,num_obs),'r--')
plot(1:num_obs,-3*sig_range*ones(1,num_obs),'r--')
title(sprintf('Range RMS = %.4e m',output.range_RMS))
ylabel('m')
subplot(2,1,2)
plot(1:num_obs, pfr_store(2,:),'.','LineWidth',1)
hold on
plot(1:num_obs,3*sig_rangerate*ones(1,num_obs),'r--')
plot(1:num_obs,-3*sig_rangerate*ones(1,num_obs),'r--')
title(sprintf('Range-Rate RMS = %.4e m/s',output.rangerate_RMS))
ylabel('m/s'),xlabel('Observation')

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
[y_sound,Fs,NBITS] = wavread('C:\Users\John\Documents\StarCraft_Sound_Pack\Protoss\Units\Artanis\patpss00');
sound(y_sound,Fs,NBITS);
end

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