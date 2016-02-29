%% SRIF

function output = SRIF(state_ap, meas_store, fo)
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
num_state_store = 6;

% Set up storage
x_est_store = zeros(num_state_store,num_obs);
prefit_range_store = zeros(num_obs,1);
state_store = zeros(num_state_store,num_obs);
state_ap_store = zeros(num_state_store,num_obs);
num_variance_store = 3;
EKF_cov_store = zeros(num_variance_store,num_obs);


STM_accum = eye(consts.state_len);
pfr_store = zeros(2,num_obs);

mystr = '';
fprintf('Observation ');

% Process observations
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
        if fo.use_DMC
            STM_obs2obs = compute_DMC_STM(fo, dt, STM_obs2obs(1:fo.important_block(1),1:fo.important_block(2)));
        end
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
        EKF_cov_store(:,ii) = diag(P(1:3,1:3));
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