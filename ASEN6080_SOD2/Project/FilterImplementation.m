%% Attitude Filter Implementations by John Clouse
%% quaternion EKF, no process noise
project_filter_params;
filter_opts.P0 = diag([0.001,0.001,0.001,0.001, 1e-5,1e-5,1e-5]);
state_ap = [1;0;0;0;0;-orbit.n;0];
EKFoutput = AttKalmanFilter(state_ap, rate_meas_data(1:200,:), filter_opts);
%% Euler EKF
project_filter_params;
filter_opts.P0 = diag([5*pi/180,5*pi/180,5*pi/180, 1e-5,1e-5,1e-5]);
state_ap = [0;0;0;0;-orbit.n;0];
% state_ap = [0;0;0;0;0;0];
state_ap = [[0;5;-3]*pi/180;0;-orbit.n;0];
filter_opts.state_length = 6;
filter_opts.important_block = [6 6];
filter_opts.propagator_opts.att_prop_opts.fcn = @rigid_body_Euler_estimation;
% EulerEKFoutput = AttKalmanFilter(state_ap, rate_meas_data(1:200,:), filter_opts);
EulerEKFoutput = AttKalmanFilter(state_ap, rate_meas_data, filter_opts);

figure; 
for ii = 1:3
subplot(3,1,ii)
plot(EulerEKFoutput.state_store(ii,:)*180/pi)
end

%% quaternion CKF, no process noise
%% Multiplicitive EKF, no process noise
%% Process noise
%% Bad ephemeris
% If the ephemeris is inaccurate for some reason, could throw things off.