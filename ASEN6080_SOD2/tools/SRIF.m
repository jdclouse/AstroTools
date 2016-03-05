%% SRIF

function output = SRIF(X0_ap, P0, meas_store, fo)
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

sig_range = 0.01; % m
sig_rangerate = 0.001; %m/s
% W = [1/(sig_range*sig_range) 0; 0 1/(sig_rangerate*sig_rangerate)];
meas_noise_cov = [(sig_range*sig_range) 0; 0 (sig_rangerate*sig_rangerate)];
V = chol(meas_noise_cov);

% init
state = X0_ap;
x_est = zeros(consts.state_len,1);
R_bar = inv(chol(P0,'upper'));
Rj = R_bar;
obs_time_last = ObsData(1,obs_t_idx);

num_state_store = 7;

% Set up storage
x_est_store = zeros(num_state_store,num_obs);
prefit_range_store = zeros(num_obs,1);
state_store = zeros(num_state_store,num_obs);
num_variance_store = 7;
cov_store = zeros(num_variance_store,num_obs);

STM_accum = eye(consts.state_len);
pfr_store = zeros(2,num_obs);

mystr = '';
fprintf('Observation ');

% Process observations
for ii = 1:num_obs
    fprintf(repmat('\b',size(mystr)))
    mystr = sprintf('%d', ii);
    fprintf(mystr);

    site_num = ObsData(ii, obs_site_idx);
    
    obs_time = ObsData(ii,obs_t_idx); 
    
    % Propagate
    % Get the state on the reference trajectory at this obs time.
    % Get STM from last obs to this one.
    dt = obs_time - obs_time_last;
    if fo.integrate_ref_state
    if ii == 1 || dt == 0
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
        
        [~,X] = ode45(@two_body_state_dot, times, last_state, ...
            fo.ode_opts, propagator_opts);
        STM_obs2obs(1:fo.important_block(1),1:fo.important_block(2)) = ...
            reshape(X(end,consts.state_len+1:end), ...
            fo.important_block(1), fo.important_block(2));

    end
    ref_state_at_obs = X(end,1:consts.state_len)';
    else
        ref_state_at_obs = fo.ref_state(ii,1:consts.state_len)';
        STM_accum_last = STM_accum;
        STM_accum = reshape(fo.ref_state(ii,consts.state_len+1:end),...
            fo.important_block(1), fo.important_block(2));
        STM_obs2obs = STM_accum_last/STM_accum;
    end
        
    % Calculate measurement deviation y
    t_obs = ObsData(ii,obs_t_idx);
    
    r_comp = compute_range_ECFsite(ref_state_at_obs(1:3),...
        site(site_num).r,theta_dot*t_obs);
    rr_comp = compute_range_rate_ECFsite(ref_state_at_obs(1:6),...
        site(site_num).r,theta_dot*t_obs, theta_dot);
    
    y1(ii) = (ObsData(ii,obs_r_idx)-r_comp);
    y2(ii) = (ObsData(ii,obs_rr_idx)-rr_comp);
    
    % Time update
    STM_accum = STM_obs2obs*STM_accum;
    x_ap = STM_obs2obs*x_est;
%     R_bar = householder( R/STM_obs2obs, 6, 6);
    R_bar = Rj/STM_obs2obs;
    b_bar = R_bar*x_ap;
    
    % H~
    consts.t = obs_time;
    consts.site = site_num;
    consts.site_r = site(site_num).r;
    H_tilda = fo.H_tilda_handle(ref_state_at_obs, consts);
%     H = H_tilda*STM_accum;
    H = H_tilda;
    
    % Measurement Update
    y = [y1(ii);y2(ii)];
    
%     xformed = householder([R_bar b_bar; V\H  V\y],8,7);
    xformed = householder([R_bar b_bar; V\H  V\y],9,8);
    Rj = xformed(1:consts.state_len,1:consts.state_len);
    bj = xformed(1:consts.state_len,end);
    
    x_est = backsub( Rj, bj, consts.state_len);
    
    % Track the last time
    obs_time_last = obs_time;    
    
    % Set up ref state for next pass
    state = ref_state_at_obs;

    pfr = y-H_tilda*x_est;
    pfr_store(:,ii) = pfr;
    
    % Store
    x_est_store(:,ii) = x_est;
    prefit_range_store(ii) = y1(ii);
    state_store(:,ii) = ...
        state(1:num_state_store) + x_est(1:num_state_store);
    Rinv = inv(Rj);
    P = Rinv*Rinv';
    cov_store(:,ii) = ...
        diag(P(1:num_variance_store,1:num_variance_store));
    
end
fprintf('\n');

RMS_accum = 0;
for ii = 1:num_obs
    RMS_accum = RMS_accum + pfr_store(:,ii)'/meas_noise_cov*pfr_store(:,ii);
end
output.RMS = sqrt(RMS_accum/num_obs);

output.rangerate_RMS = sqrt(sum(pfr_store(2,:).*pfr_store(2,:))/num_obs);
output.range_RMS = sqrt(sum(pfr_store(1,:).*pfr_store(1,:))/num_obs);
output.pfr_store = pfr_store;
output.prefit_range_store = prefit_range_store;
output.cov_store = cov_store;
output.state_store = state_store;
output.x_est_store = x_est_store;
% output.state_ap_store = state_ap_store;