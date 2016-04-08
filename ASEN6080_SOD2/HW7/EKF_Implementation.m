filter_opts.use_EKF = 1;
sig_range = 0.005; % km
sig_rangerate = 0.5*1e-6; %m/s
% W = [1/(sig_range*sig_range) 0; 0 1/(sig_rangerate*sig_rangerate)];
filter_opts.R = [(sig_range*sig_range) 0; 0 (sig_rangerate*sig_rangerate)];
filter_opts.P0 = P;

EKFoutput = KalmanFilter(state_ap, ObsMassaged, filter_opts);