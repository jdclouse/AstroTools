function output = UnscentedKalmanFilter(state_ap, P_ap, meas_store, fo)

ObsData = meas_store;
[num_obs, ~] = size(ObsData);
% num_obs = 500;
obs_r_idx = 3;
obs_rr_idx = 4;
obs_t_idx = 2;
obs_site_idx = 1;

site = fo.site;
theta_dot = fo.theta_dot;
fo.ode_opts = odeset('RelTol', 1e-10, 'AbsTol', 1e-10);

L = fo.UKF.L; 
alpha = fo.UKF.alpha; 
beta = fo.UKF.beta; 

K = 3-L; %FIXME options
lambda = alpha*alpha*(L+K)-L;
w_0_m = lambda/(L+lambda);
w_0_c = w_0_m + 1-alpha*alpha+beta;
w_i = 1/(2*(L+lambda));

num_sig_pts = (2*L+1);

% weights in matrix to get rid of looping.
mean_weights = repmat([w_0_m repmat(w_i,1,num_sig_pts-1)],L,1);

P = P_ap;
X_est = state_ap;
% R;
% FIXME
sig_range = 0.01; % m
sig_rangerate = 0.001; %m/s
R = [(sig_range*sig_range) 0; 0 (sig_rangerate*sig_rangerate)];
last_obs_time = 0;

gamma = alpha*sqrt(3);

% Storage
X_est_store = zeros(L,num_obs);
Pt_store = zeros(L,num_obs);
pfr_store = zeros(2,num_obs);

mystr = '';
fprintf('Observation ');

% Run the UKF
for ii = 1:num_obs
    % "status bar"
    fprintf(repmat('\b',size(mystr)))
    mystr = sprintf('%d', ii);
    fprintf(mystr);
    
    % Get everything needed from observation data
    t_obs = ObsData(ii,obs_t_idx);
    site_num = ObsData(ii, obs_site_idx);
    y_obs = [ObsData(ii,obs_r_idx);...
        ObsData(ii,obs_rr_idx)];

    % Set up sig pts
    sqrtP = sqrtm(P);
    sig_pts = [X_est, repmat(X_est,1,L) + gamma*sqrtP, ...
        repmat(X_est,1,L) - gamma*sqrtP];
    
    dt = t_obs - last_obs_time;
    if dt > 0

        % Propagate
        prop_state_last = reshape(sig_pts, L*num_sig_pts, 1);
        times = [0 dt];

        [~,X] = ode45(@UKF_two_body_state_dot, times, prop_state_last, ...
            fo.ode_opts, fo.propagator_opts);

        sig_pts_new = reshape(X(end,:)',L,num_sig_pts);
    else
        
        sig_pts_new = sig_pts;
%         X_ap = X_est;
%         Pbar_t = P;
    end
    
    % Time update
    X_ap = w_0_m*sig_pts_new(:,1);
    for jj = 2:num_sig_pts
        X_ap = X_ap + w_i*sig_pts_new(:,jj);
    end
    X_ap = sum(mean_weights.*sig_pts_new,2);

    Pbar_t = w_0_c*(sig_pts_new(:,1) - X_ap)*(sig_pts_new(:,1) - X_ap)';
    for jj = 2:num_sig_pts
        Pbar_t = Pbar_t + ...
            w_i*(sig_pts_new(:,jj) - X_ap)*(sig_pts_new(:,jj) - X_ap)';
    end

    if fo.use_SNC ...
            && dt < fo.SNC_meas_separation_threshold

        % Compute Q. Use the state from last obs for any transformations.
        Q = compute_SNC_Q(fo, X_est(1:L));
        Gamma = fo.SNC_Gamma(dt);

        Pbar_t = Pbar_t + Gamma*Q*Gamma'; % Add the process noise
    end

    % Recompute sig-pts
    % These are basically a priori sig pts
    sqrt_Pbar = sqrtm(Pbar_t);
    sig_pts = [X_ap, repmat(X_ap,1,L) + gamma*sqrt_Pbar, ...
        repmat(X_ap,1,L) - gamma*sqrt_Pbar];
    
    % Computed measurements and weighted averages
    y_mean = zeros(2,1);
    y_sig_comp = zeros(2,2*L+1);
    for jj = 1:num_sig_pts
        r_comp = compute_range_ECFsite(sig_pts(1:3,jj),...
            site(site_num).r*1e3,theta_dot*t_obs);
        rr_comp = compute_range_rate_ECFsite(sig_pts(1:6,jj),...
            site(site_num).r*1e3,theta_dot*t_obs, theta_dot);
        y_sig_comp(:,jj) = [r_comp; rr_comp];
        if jj == 1
            y_mean = y_mean + w_0_m*y_sig_comp(:,jj);
        else
            y_mean = y_mean + w_i*y_sig_comp(:,jj);
        end
    end
    
    % Innovation and cross-correlation
    Pyy = R + w_0_c*(y_sig_comp(:,1)-y_mean)*(y_sig_comp(:,1)-y_mean)';
    Pxy = w_0_c*(sig_pts_new(:,1) - X_ap)*(y_sig_comp(:,1)-y_mean)';
    for jj = 2:num_sig_pts
        Pyy = Pyy + w_i*(y_sig_comp(:,jj)-y_mean)*(y_sig_comp(:,jj)-y_mean)';
        Pxy = Pxy + ...
            w_i*(sig_pts_new(:,jj) - X_ap)*(y_sig_comp(:,jj)-y_mean)';
    end
    
    % Measurement update
    Kt = Pxy/Pyy;
    X_est = X_ap + Kt*(y_obs - y_mean);
%     P = Pbar_t - Kt*Pyy*Kt';
    P = Pbar_t - Pxy*Kt'; % this removes an inversion
    
    % Post-fit residuals
%     if alpha == 1
    y_mean = zeros(2,1);
    sqrtP = sqrtm(P);
    sig_pts = [X_est, repmat(X_est,1,L) + gamma*sqrtP, ...
        repmat(X_est,1,L) - gamma*sqrtP];
    for jj = 1:num_sig_pts
        r_comp = compute_range_ECFsite(sig_pts(1:3,jj),...
            site(site_num).r*1e3,theta_dot*t_obs);
        rr_comp = compute_range_rate_ECFsite(sig_pts(1:6,jj),...
            site(site_num).r*1e3,theta_dot*t_obs, theta_dot);
        y_sig_comp = [r_comp; rr_comp];
        if jj == 1
            y_mean = y_mean + w_0_m*y_sig_comp;
        else
            y_mean = y_mean + w_i*y_sig_comp;
        end
    end
%     end
    pfr = y_obs - y_mean;
    if ~isreal(pfr)
        zzzzz = 1;
    end
    
    % Reset
    last_obs_time = t_obs;
    
    % Store everything
    X_est_store(:,ii) = X_est;
    Pt_store(:,ii) = diag(P);
    pfr_store(:,ii) = pfr;
end

% set the output
output.X_est_store = X_est_store;
output.Pt_store = Pt_store;
output.pfr_store = pfr_store;
output.range_RMS = sqrt(sum(pfr_store(1,:).*pfr_store(1,:))/num_obs);
output.rangerate_RMS = sqrt(sum(pfr_store(2,:).*pfr_store(2,:))/num_obs);

end