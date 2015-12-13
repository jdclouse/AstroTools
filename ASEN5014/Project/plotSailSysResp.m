function  plotSailSysResp( analysis_set,y,t,K,r,Ts,torque_tmax )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
% Coning angle plot

figWidth = 1120; % pixels
figHeight = 840; % pixels
r2d = 180/pi;
gimbal_angle_lim = pi/6*r2d;

ctrl1Obs_Alpha_plot = figure('Position', [0, 0, figWidth, figHeight]);
plot(t/3600,y(:,1)*r2d)
ylabel('\alpha (deg)')
xlabel('Time (hr)')
title(...
    [sprintf('Sun Angle, Step Reference at %.1f degrees', r*180/pi)...
    ', ' analysis_set]);
hold on
plot([Ts Ts]/3600,[0, r*1.10]*r2d,'r')
plot([t(1),t(end)]/3600,[1 1]*r*1.1*r2d,'r-.')
plot([t(1),t(end)]/3600,[1 1]*r*1.05*r2d,'r--')
plot([t(1),t(end)]/3600,[1 1]*r*.95*r2d,'r--')
legend('\alpha','T_{settle}','10% Overshoot', '5% Settling',...
    'Location','SouthEast')
print(['Report/' analysis_set '_Alpha'],'-dpng')

% Control torque
u = [];
for ii = 1:length(t)
    u(ii) = r-K*y(ii,:)';
end
ctrl1Obs_Torque_plot = figure('Position', [0, 0, figWidth, figHeight]);
% torque_tmax = 600;
plot(t(t<torque_tmax),u(t<torque_tmax))
title([sprintf('Gimbal Torque, %.0f-degree Step',r*180/pi)...
    ', ' analysis_set])
ylabel('T_g (N-m)')
xlabel('Time (s)')
print(['Report/' analysis_set '_Torque'],'-dpng')

%Gimbal angle plot
ctrl1Obs_GimbalAng_plot = figure('Position', [0, 0, figWidth, figHeight]);
plot(t/3600,y(:,3)*r2d)
hold on
plot([t(1) t(end)]/3600,[1 1]*gimbal_angle_lim, 'r--')
plot([t(1) t(end)]/3600,[1 1]*-gimbal_angle_lim, 'r--')
title([sprintf('Gimbal Angle, %.0f-degree Step',r*180/pi)...
    ', ' analysis_set])
ylabel('\delta (deg)')
ylim([-gimbal_angle_lim-1, gimbal_angle_lim+1])
xlabel('Time (hr)')
legend('\delta','Actuator Limit',...
    'Location','SouthEast')
print(['Report/' analysis_set '_Delta'],'-dpng')

end

