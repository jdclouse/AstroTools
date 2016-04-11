filter_opts.use_EKF = 1;
filter_opts.use_SNC = 1;
filter_opts.SNC_Q = eye(3)*1e-17;
sig_range = 0.005; % km
sig_rangerate = 0.5*1e-6; %m/s
% W = [1/(sig_range*sig_range) 0; 0 1/(sig_rangerate*sig_rangerate)];
filter_opts.R = [(sig_range*sig_range) 0; 0 (sig_rangerate*sig_rangerate)];
filter_opts.P0 = P;
filter_opts.use_smoother = 1;
% (1:1000,:)
tic
EKFoutput = KalmanFilter(state_ap, ObsMassaged(1:12000,:), filter_opts);
toc

if filter_opts.use_SNC
    plot_title = sprintf('EKF w/ SNC, Q=%.1e', filter_opts.SNC_Q(1,1));
else
    plot_title = 'EKF w/o SNC';
end
residual_plot(EKFoutput.pfr_store, [0.005, 0.5*1e-6], plot_title)
residual_plot(EKFoutput.pfr_store, [0.005, 0.5*1e-6], plot_title, ObsMassaged(:,2))
EKFoutput.smoothed_state_store(:,1)

ell_plot = figure;
plot_handles1 = [];
output_state = EKFoutput.state_store(:,end);
final_P = EKFoutput.final_P;
color_num = 1;
BPlaneTarget_HW7;
output_state = EKFoutputorig.state_store(:,end);
final_P = EKFoutputorig.final_P;
color_num = 2;
BPlaneTarget_HW7;

figure(ell_plot);
legend(plot_handles1, 'Full obs set');
xlabel('T (km)')
ylabel('R (km)')
title('3\sigma covariance ellipse on B-plane')
set(gca,'YDir','reverse')