% UKF
filter_opts.propagator_opts.UKF.L = 7;
filter_opts.UKF.L = 7; 
filter_opts.UKF.alpha = 1; 
filter_opts.UKF.beta = 2; 
sig_range = 0.005; % km
sig_rangerate = 0.5*1e-6; %m/s
filter_opts.R = [(sig_range*sig_range) 0; 0 (sig_rangerate*sig_rangerate)];

filter_opts.use_SNC = 1;
filter_opts.SNC_Q = eye(3)*1e-5;
filter_opts.SNC_Gamma = @(dt) [dt*dt/2 0 0;...
            0 dt*dt/2 0;...
            0 0 dt*dt/2;...
            dt 0 0;...
            0 dt 0;...
            0 0 dt;
            0 0 0];
filter_opts.propagator_opts.J3.use = 0;
tic
UKFstorage = UnscentedKalmanFilter(state_ap, P, ObsMassaged, filter_opts);
toc

residual_plot(UKFstorage.pfr_store, [0.005, 0.5*1e-6], '1 iter')
