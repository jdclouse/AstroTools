filter_opts.use_EKF = 1;
filter_opts.use_SNC = 1;
filter_opts.SNC_Q = eye(3)*1e-13;
sig_range = 0.005; % km
sig_rangerate = 0.5*1e-6; %m/s
% W = [1/(sig_range*sig_range) 0; 0 1/(sig_rangerate*sig_rangerate)];
filter_opts.R = [(sig_range*sig_range) 0; 0 (sig_rangerate*sig_rangerate)];
filter_opts.P0 = P;
% (1:1000,:)
tic
EKFoutput = KalmanFilter(state_ap, ObsMassaged, filter_opts);
toc

if filter_opts.use_SNC
    plot_title = sprintf('EKF w/ SNC, Q=%.1e', filter_opts.SNC_Q(1,1));
else
    plot_title = 'EKF w/o SNC';
end
residual_plot(EKFoutput.pfr_store, [0.005, 0.5*1e-6], plot_title)
residual_plot(EKFoutput.pfr_store, [0.005, 0.5*1e-6], plot_title, ObsMassaged(:,2))