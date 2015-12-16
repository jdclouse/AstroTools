function  plotSailSysResp( analysis_set,y,t,K,r,Ts,torque_tmax, use_title )
%plotSailSysResp All the controller plots


y_size = size(y);
if length(y_size) == 2
    num_y = 1;
elseif length(y_size) == 3
    num_y = y_size(end);
end

figWidth = 1120; % pixels
figHeight = 840; % pixels
r2d = 180/pi;
gimbal_angle_lim = pi/6*r2d;
lw = 1;
fs = 20;

% Coning angle plot
ctrl1Obs_Alpha_plot = figure('Position', [0, 0, figWidth, figHeight]);

for jj = 1:num_y
h_plot = plot(t/3600,y(:,1,jj)*r2d,'LineWidth',lw);
hold on
end
ylabel('\alpha (deg)', 'FontSize', fs)
xlabel('Time (hr)', 'FontSize', fs)
if use_title
title(...
    [sprintf('Sun Angle, Step Reference at %.1f degrees', r*180/pi)...
    ', ' analysis_set], 'FontSize', fs);
end
hold on
h_Ts = plot([Ts Ts]/3600,[0, r*1.10]*r2d,'r','LineWidth',lw);
h_OS = plot([t(1),t(end)]/3600,[1 1]*r*1.1*r2d,'r-.','LineWidth',lw);
h_hi = plot([t(1),t(end)]/3600,[1 1]*r*1.05*r2d,'r--','LineWidth',lw);
h_lo = plot([t(1),t(end)]/3600,[1 1]*r*.95*r2d,'r--','LineWidth',lw);
set(gca, 'FontSize', fs)
pl = legend([h_plot h_Ts h_OS h_hi],...
    '\alpha','T_{settle}','10% Overshoot', '5% Settling',...
    'Location','SouthEast');
set(pl, 'FontSize', fs)
print(['Report/' analysis_set '_Alpha'],'-dpng')

% Control torque
ctrl1Obs_Torque_plot = figure('Position', [0, 0, figWidth, figHeight]);
for jj = 1:num_y
u = [];
for ii = 1:length(t)
    u(ii) = r-K*y(ii,:,jj)';
end
plot(t(t<torque_tmax),u(t<torque_tmax),'LineWidth',lw)
hold on
end
if use_title
title([sprintf('Gimbal Torque, %.0f-degree Step',r*180/pi)...
    ', ' analysis_set], 'FontSize', fs)
end
set(gca, 'FontSize', fs)
ylabel('T_g (N-m)', 'FontSize', fs)
xlabel('Time (s)', 'FontSize', fs)
print(['Report/' analysis_set '_Torque'],'-dpng')

%Gimbal angle plot
ctrl1Obs_GimbalAng_plot = figure('Position', [0, 0, figWidth, figHeight]);
for jj = 1:num_y
h_plot = plot(t/3600,y(:,3,jj)*r2d,'LineWidth',lw);
hold on
end
h_hi = plot([t(1) t(end)]/3600,[1 1]*gimbal_angle_lim, 'r--','LineWidth',lw);
h_lo = plot([t(1) t(end)]/3600,[1 1]*-gimbal_angle_lim, 'r--','LineWidth',lw);
if use_title
title([sprintf('Gimbal Angle, %.0f-degree Step',r*180/pi)...
    ', ' analysis_set], 'FontSize', fs)
end
set(gca, 'FontSize', fs)
ylabel('\delta (deg)', 'FontSize', fs)
ylim([-gimbal_angle_lim-1, gimbal_angle_lim+1])
xlabel('Time (hr)', 'FontSize', fs)
pl = legend([h_plot, h_hi], '\delta','Actuator Limit',...
    'Location','SouthEast');
set(pl, 'FontSize', fs)
print(['Report/' analysis_set '_Delta'],'-dpng')

end

