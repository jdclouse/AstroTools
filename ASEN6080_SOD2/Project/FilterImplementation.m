%% Attitude Filter Implementations by John Clouse
%% quaternion EKF, no process noise
project_filter_params;
filter_opts.P0 = diag([0.001,0.001,0.001,0.001, 1e-5,1e-5,1e-5]);
state_ap = [1;0;0;0;0;-orbit.n;0];
EKFoutput = AttKalmanFilter(state_ap, rate_meas_data(1:200,:), filter_opts);
%% Euler CKF
project_filter_params;
filter_opts.P0 = diag([(5*pi/180)^2,(5*pi/180)^2,(5*pi/180)^2, 1e-5,1e-5,1e-5]);
% state_ap = [0;0;0;0;-orbit.n;0];
% state_ap = [0;0;0;0;0;0];
state_ap = [[0;5;-3]*pi/180;0;-orbit.n;0];
% state_ap = [[0;4;-2]*pi/180;0;-orbit.n;0];
filter_opts.state_length = 6;
filter_opts.important_block = [6 6];
filter_opts.propagator_opts.att_prop_opts.fcn = @rigid_body_Euler_estimation;

filter_opts.use_EKF = 0;
filter_opts.use_SNC = 1;
% 1e-10 fit the data way too much, good example though
filter_opts.SNC_Q = eye(3)*1e-14;
% EulerCKFoutput = AttKalmanFilter(state_ap, rate_meas_data(1:500,:), filter_opts);
EulerCKFoutput = AttKalmanFilter(state_ap, rate_meas_data, filter_opts);
% 
% figure; 
% for ii = 1:3
% subplot(3,1,ii)
% plot(EulerCKFoutput.state_store(ii,:)*180/pi)
% end
% subplot(3,1,1);
% title(sprintf('CKF Angles, SNC = %.1e',filter_opts.SNC_Q(1)))
% % rates
% figure; 
% for ii = 1:3
% subplot(3,1,ii)
% plot(EulerCKFoutput.state_store(ii+3,:)*180/pi)
% hold on
% plot(rate_meas_data(:,1+ii)*180/pi, 'r')
% end
% subplot(3,1,1);
% title(sprintf('CKF Rates, SNC = %.1e',filter_opts.SNC_Q(1)))
% 
% % Error
% figure; 
% for ii = 1:3
% subplot(3,1,ii)
% plot((EulerCKFoutput.state_store(ii,:)-X_out_angles(:,ii+6)')*180/pi)
% end
% subplot(3,1,1);
% title(sprintf('CKF Angle Error, SNC = %.1e',filter_opts.SNC_Q(1)))

filter_output = EulerCKFoutput;
case_title = 'CKF - Rates only';
resid_plot_units = {'$\delta\dot{\psi}$ (rad/s)', '$\delta\dot{\theta}$ (rad/s)', '$\delta\dot{\phi}$ (rad/s)'};
plot_results;

%% Euler CKF -- with yaw
project_filter_params;
filter_opts.P0 = diag([5*pi/180,5*pi/180,5*pi/180, 1e-5,1e-5,1e-5]);
state_ap = [0;0;0;0;-orbit.n;0];
% state_ap = [[0;5;-3]*pi/180;0;-orbit.n;0];
% state_ap = [[0;4;-2]*pi/180;0;-orbit.n;0];
filter_opts.state_length = 6;
filter_opts.important_block = [6 6];
filter_opts.propagator_opts.att_prop_opts.fcn = @rigid_body_Euler_estimation;

filter_opts.R = eye(4)*(1e-6)^2;
filter_opts.R(4,4) = (0.5/3*pi/180)^2;

filter_opts.use_EKF = 0;
filter_opts.use_SNC = 1;
% 1e-10 fit the data way too much, good example though
filter_opts.SNC_Q = eye(3)*1e-14;
% EulerCKFoutput_yaw = AttKalmanFilter(state_ap, yaw_and_rate_meas_data(1:500,:), filter_opts);
EulerCKFoutput_yaw = AttKalmanFilter(state_ap, yaw_and_rate_meas_data, filter_opts);

figure; 
for ii = 1:3
subplot(3,1,ii)
plot(EulerCKFoutput_yaw.state_store(ii,:)*180/pi)
hold on
plot(X_out_angles(:,ii+6)'*180/pi, 'r')
end
subplot(3,1,1);
title(sprintf('CKF Angles, SNC = %.1e',filter_opts.SNC_Q(1)))
% rates
figure; 
for ii = 1:3
subplot(3,1,ii)
plot(EulerCKFoutput_yaw.state_store(ii+3,:)*180/pi)
hold on
plot(rate_meas_data(:,1+ii)*180/pi, 'r')
end
subplot(3,1,1);
title(sprintf('CKF Rates, SNC = %.1e',filter_opts.SNC_Q(1)))

% Error
figure; 
for ii = 1:3
subplot(3,1,ii)
plot((EulerCKFoutput_yaw.state_store(ii,:)-X_out_angles(:,ii+6)')*180/pi)
end
subplot(3,1,1);
title(sprintf('CKF Angle Error, SNC = %.1e',filter_opts.SNC_Q(1)))

plot_att_resid(EulerCKFoutput_yaw.pfr_store, filter_opts.R, {'rad/s', 'rad/s', 'rad/s', 'rad'})

%% Euler EKF -- yaw
project_filter_params;
filter_opts.P0 = diag([(5*pi/180)^2,(5*pi/180)^2,(5*pi/180)^2, 1e-5,1e-5,1e-5]);
state_ap = [0;0;0;0;-orbit.n;0];
% state_ap = [0;0;0;0;0;0];
% state_ap = [[0;5;-3]*pi/180;0;-orbit.n;0];
state_ap = [[0;4;-2]*pi/180;0;-orbit.n;0];
filter_opts.state_length = 6;
filter_opts.important_block = [6 6];
filter_opts.propagator_opts.att_prop_opts.fcn = @rigid_body_Euler_estimation;

filter_opts.R = eye(4)*(1e-6)^2;
filter_opts.R(4,4) = (0.5/3*pi/180)^2;

filter_opts.use_EKF = 1;
filter_opts.EKF_switchover = 1000;
filter_opts.use_SNC = 1;
% 1e-10 fit the data way too much, good example though
filter_opts.SNC_Q = eye(3)*1e-14;
% EulerEKFoutput = AttKalmanFilter(state_ap, yaw_and_rate_meas_data(1:2000,:), filter_opts);
EulerEKFoutput = AttKalmanFilter(state_ap, yaw_and_rate_meas_data, filter_opts);

filter_output = EulerEKFoutput;
case_title = 'EKF - Yaw and Rate Meas';
resid_plot_units = {'$\delta\dot{\psi}$ (rad/s)', '$\delta\dot{\theta}$ (rad/s)', '$\delta\dot{\phi}$ (rad/s)', '$\psi$ (rad)'};
plot_results;

%% Euler EKF -- YPR measurements
project_filter_params;
filter_opts.P0 = diag([5*pi/180,5*pi/180,5*pi/180, 1e-5,1e-5,1e-5]);
state_ap = [0;0;0;0;-orbit.n;0];
% state_ap = [0;0;0;0;0;0];
% state_ap = [[0;5;-3]*pi/180;0;-orbit.n;0];
state_ap = [[0;4;-2]*pi/180;0;-orbit.n;0];
filter_opts.state_length = 6;
filter_opts.important_block = [6 6];
filter_opts.propagator_opts.att_prop_opts.fcn = @rigid_body_Euler_estimation;

filter_opts.R = eye(3)*(0.5/3*pi/180)^2;

filter_opts.use_EKF = 1;
filter_opts.EKF_switchover = 1000;
filter_opts.use_SNC = 1;
% 1e-13 had residuals that were too big.
filter_opts.SNC_Q = eye(3)*1e-16;

filter_opts.YPR = true;
EulerEKFoutput_YPR = AttKalmanFilter(state_ap, YPR_meas_data(1:2000,:), filter_opts);
% EulerEKFoutput_YPR = AttKalmanFilter(state_ap, YPR_meas_data, filter_opts);

filter_output = EulerEKFoutput_YPR;
case_title = 'EKF - YPR';
resid_plot_units = {'rad', 'rad', 'rad'};
plot_results;

%% Euler EKF -- YPR and rate measurements
project_filter_params;
filter_opts.P0 = diag([5*pi/180,5*pi/180,5*pi/180, 1e-5,1e-5,1e-5]);
state_ap = [0;0;0;0;-orbit.n;0];
% state_ap = [0;0;0;0;0;0];
% state_ap = [[0;5;-3]*pi/180;0;-orbit.n;0];
% state_ap = [[0;4;-2]*pi/180;0;-orbit.n;0];
filter_opts.state_length = 6;
filter_opts.important_block = [6 6];
filter_opts.propagator_opts.att_prop_opts.fcn = @rigid_body_Euler_estimation;

filter_opts.R = eye(6)*(0.5/3*pi/180)^2;
filter_opts.R(4:6,4:6) = eye(3)*(1e-6)^2;

filter_opts.use_EKF = 1;
filter_opts.EKF_switchover = 200;
filter_opts.use_SNC = 1;
% 1e-13 had residuals that were too big.
filter_opts.SNC_Q = eye(3)*1e-14;

filter_opts.YPR = true;
filter_opts.YPR_rates = true;
% EulerEKFoutput_YPR_rates = AttKalmanFilter(state_ap, YPR_meas_data_and_rates(1:2000,:), filter_opts);
EulerEKFoutput_YPR_rates = AttKalmanFilter(state_ap, YPR_meas_data_and_rates, filter_opts);

filter_output = EulerEKFoutput_YPR_rates;
case_title = 'EKF - YPR';
resid_plot_units = {'rad', 'rad', 'rad', 'rad/s', 'rad/s', 'rad/s'};
plot_results;

%% Euler EKF -- YPR and rate measurements, SNC Gamma changed
project_filter_params;
filter_opts.P0 = diag([5*pi/180,5*pi/180,5*pi/180, 1e-5,1e-5,1e-5]);
state_ap = [0;0;0;0;-orbit.n;0];
% state_ap = [0;0;0;0;0;0];
% state_ap = [[0;5;-3]*pi/180;0;-orbit.n;0];
% state_ap = [[0;4;-2]*pi/180;0;-orbit.n;0];
filter_opts.state_length = 6;
filter_opts.important_block = [6 6];
filter_opts.propagator_opts.att_prop_opts.fcn = @rigid_body_Euler_estimation;

filter_opts.R = eye(6)*(0.5/3*pi/180)^2;
filter_opts.R(4:6,4:6) = eye(3)*(1e-6)^2;

filter_opts.use_EKF = 1;
filter_opts.EKF_switchover = 200;
filter_opts.use_SNC = 1;
filter_opts.SNC_Q = eye(3)*1e-13;
filter_opts.SNC_Gamma = @(dt) [dt 0 0;...
            0 dt 0;...
            0 0 dt;...
            1 0 0;...
            0 1 0;...
            0 0 1];

filter_opts.YPR = true;
filter_opts.YPR_rates = true;
% EulerEKFoutput_YPR_rates_GammaChange = AttKalmanFilter(state_ap, YPR_meas_data_and_rates(1:2000,:), filter_opts);
EulerEKFoutput_YPR_rates_GammaChange = AttKalmanFilter(state_ap, YPR_meas_data_and_rates, filter_opts);

filter_output = EulerEKFoutput_YPR_rates_GammaChange;
case_title = 'EKF - YPR';
resid_plot_units = {'rad', 'rad', 'rad', 'rad/s', 'rad/s', 'rad/s'};
plot_results;
%% Euler EKF -- YPR and rate measurements, wrong altitude
project_filter_params;

n_off = sqrt(Earth.mu/(orbit.a+10)^3);

[r_init_off, v_init_off] = OE2cart(orbit.a+10, orbit.e, orbit.i, orbit.RAAN, ...
    orbit.w, orbit.f, Earth.mu);
filter_opts.PV_state = [r_init_off; v_init_off];

filter_opts.P0 = diag([5*pi/180,5*pi/180,5*pi/180, 1e-5,1e-5,1e-5]);
state_ap = [0;0;0;0;-n_off;0];
% state_ap = [0;0;0;0;0;0];
% state_ap = [[0;5;-3]*pi/180;0;-orbit.n;0];
% state_ap = [[0;4;-2]*pi/180;0;-orbit.n;0];
filter_opts.state_length = 6;
filter_opts.important_block = [6 6];
filter_opts.propagator_opts.att_prop_opts.fcn = @rigid_body_Euler_estimation;

filter_opts.R = eye(6)*(0.5/3*pi/180)^2;
filter_opts.R(4:6,4:6) = eye(3)*(1e-6)^2;

filter_opts.use_EKF = 1;
filter_opts.EKF_switchover = 200;
filter_opts.use_SNC = 1;
% 1e-13 had residuals that were too big.
filter_opts.SNC_Q = eye(3)*1e-14;

filter_opts.YPR = true;
filter_opts.YPR_rates = true;

% EulerEKFoutput_YPR_rates_alt_off = AttKalmanFilter(state_ap, YPR_meas_data_and_rates(1:2000,:), filter_opts);
EulerEKFoutput_YPR_rates_alt_off = AttKalmanFilter(state_ap, YPR_meas_data_and_rates, filter_opts);

filter_output = EulerEKFoutput_YPR_rates_alt_off;
case_title = 'EKF - YPR, alt off';
resid_plot_units = {'rad', 'rad', 'rad', 'rad/s', 'rad/s', 'rad/s'};
plot_results;


%% Euler EKF
project_filter_params;
filter_opts.P0 = diag([5*pi/180,5*pi/180,5*pi/180, 1e-5,1e-5,1e-5]);
state_ap = [0;0;0;0;-orbit.n;0];
% state_ap = [0;0;0;0;0;0];
% state_ap = [[0;5;-3]*pi/180;0;-orbit.n;0];
filter_opts.state_length = 6;
filter_opts.important_block = [6 6];
filter_opts.propagator_opts.att_prop_opts.fcn = @rigid_body_Euler_estimation;

filter_opts.use_EKF = 1;
filter_opts.EKF_switchover = 1000;
filter_opts.use_SNC = 1;
% 1e-10 fit the data way too much, good example though
filter_opts.SNC_Q = eye(3)*1e-14;
% EulerEKFoutput = AttKalmanFilter(state_ap, rate_meas_data(1:500,:), filter_opts);
EulerEKFoutput = AttKalmanFilter(state_ap, rate_meas_data, filter_opts);

figure; 
for ii = 1:3
subplot(3,1,ii)
plot(EulerEKFoutput.state_store(ii,:)*180/pi)
end
subplot(3,1,1);
title(sprintf('EKF Angles, SNC = %.1e',filter_opts.SNC_Q(1)))
% rates
figure; 
for ii = 1:3
subplot(3,1,ii)
plot(EulerEKFoutput.state_store(ii+3,:)*180/pi)
hold on
plot(rate_meas_data(:,1+ii)*180/pi, 'r')
end
subplot(3,1,1);
title(sprintf('EKF Rates, SNC = %.1e',filter_opts.SNC_Q(1)))

% Error
figure; 
for ii = 1:3
subplot(3,1,ii)
plot((EulerEKFoutput.state_store(ii,:)-X_out_angles(:,ii+6)')*180/pi)
end
subplot(3,1,1);
title(sprintf('EKF Angle Error, SNC = %.1e',filter_opts.SNC_Q(1)))

%% quaternion CKF, no process noise
%% Multiplicitive EKF, no process noise
%% Process noise
%% Bad ephemeris
% If the ephemeris is inaccurate for some reason, could throw things off.