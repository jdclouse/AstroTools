%% HW8 Problem 1
fprintf('\n');
clearvars -except function_list hw_pub toolsPath 
close all
CelestialConstants; % import useful constants

r = 6378+800; %km
p = r;
e = 0;
n = sqrt(Earth.mu/r^3)*180/pi; %rad/s

i_vec = 0:0.01:180;
nodal_regression_rate = -3*n*Earth.R^2*Earth.J2/2/p^2*cosd(i_vec)*day2sec;

figure('Position',[0 0 hw_pub.figWidth hw_pub.figHeight])
plot(i_vec, nodal_regression_rate,'LineWidth',2)
xlabel('Inclination (deg)')
ylabel('$\dot{\Omega}$ (deg/day)','interpreter','latex')
title('Daily Nodal Regression')
grid on

sun_sync_node_rate = 360/365.2421897; %deg/day
sun_sync_incl = acosd(sun_sync_node_rate*2*p^2/...
    (-3*n*Earth.R^2*Earth.J2*day2sec));

hold on
plot([sun_sync_incl],[sun_sync_node_rate],'ro','LineWidth',2)

fprintf('Inclination for this orbit to be sun-synchronous: %.2f deg\n',...
    sun_sync_incl) 