function output = UnscentedKalmanFilter(state_ap, meas_store, fo)

ObsData = meas_store;
obs_r_idx = 3;
obs_rr_idx = 4;
obs_t_idx = 2;
obs_site_idx = 1;

L = 6;
alpha = 1;
beta = 2;
K = 3-L;
lambda = alpha*alpha*(L+K)-L;
w_0_m = lambda/(L+lambda);
w_0_c = w_0_m + 1-alpha*alpha+beta;
w_i = 1/(2*(L+lambda));

P;
X_est = state_ap;
R;
last_obs_time = 0;

gamma = alpha*sqrt(3);

% Storage
X_est_store = zeros(num_state_store,num_obs);
Pt_store = zeros(num_state_store,num_obs);

for ii = 1:num_obs
    
    % Get everything needed from observation data
    obs_time = ObsData(ii,obs_t_idx);
    site_num = ObsData(ii, obs_site_idx);
    y_obs = [ObsData(ii,obs_r_idx);...
        ObsData(ii,obs_rr_idx)];
    
    dt = obs_time - last_obs_time;
    if dt > 0
        % Set up sig pts
        sqrtP = sqrtm(P);
        sig_pts = [X_est, repmat(X_est,1,L) + gamma*sqrtP, ...
            repmat(X_est,1,L) - gamma*sqrtP];

        % Propagate
        prop_state_last = reshape(sig_pts, L*(2*L+1), 1);
        times = [0 dt];

        [~,X] = ode45(@two_body, times, prop_state_last, fo.ode_opts, ...
            fo.ukf_prop_opts);

        sig_pts_new = reshape(X(end,:)',L,(2*L+1));
    else
        sig_pts_new = sig_pts;
    end
    
    % Time update
    X_ap = w_0_m*sig_pts_new(1,:);
    for jj = 2:(2*L+1)
        X_ap = X_ap + w_i*sig_pts_new(jj,:);
    end
    
    Pbar_t = w_0_c*(sig_pts_new(:,1) - X_ap)*(sig_pts_new(:,1) - X_ap)';
    for jj = 2:(2*L+1)
        Pbar_t = Pbar_t + ...
            w_i*(sig_pts_new(:,jj) - X_ap)*(sig_pts_new(:,jj) - X_ap)';
    end
    
    % Recompute sig-pts
    sqrt_Pbar = sqrtm(Pbar_t);
    sig_pts = [X_ap, repmat(X_ap,1,L) + gamma*sqrt_Pbar, ...
        repmat(X_ap,1,L) - gamma*sqrt_Pbar];
    
    % Computed measurements and weighted averages
    y_mean = zeros(2,1);
    y_meas = zeros(2,2*L+1);
    for jj = 1:(2*L+1)
        r_comp = compute_range_ECFsite(sig_pts(1:3,jj),...
            site(site_num).r*1e3,theta_dot*t_obs);
        rr_comp = compute_range_rate_ECFsite(ref_state_at_obs(1:6,jj),...
            site(site_num).r*1e3,theta_dot*t_obs, theta_dot);
        y_meas(:,jj) = [r_comp; rr_comp];
        if jj == 1
            y_mean = y_mean + w_0_m*y_meas(:,jj);
        else
            y_mean = y_mean + w_i*y_meas(:,jj);
        end
    end
    
    % Innovation and cross-correlation
    Pyy = R + w_0_c*(y_meas(1,jj)-y_mean)*(y_meas(1,jj)-y_mean)';
    Pxy = w_0_c*(sig_pts_new(:,jj) - X_ap)*(y_meas(1,jj)-y_mean)';
    
    % Measurement update
    Kt = Pxy/Pyy;
    X_est = X_ap + Kt*(y_obs - y_mean);
    Pt = Pbar_t - Kt*Pyy*Kt';
    
    % Store everything
    X_est_store(:,ii) = X_est;
    Pt_store(:,ii) = diag(Pt);
end

output.X_est_store = X_est_store;
output.Pt_store = Pt_store;

end