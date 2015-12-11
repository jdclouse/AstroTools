%% ASEN 5050 Project: Solar Sail Trajectories
%% Setup
clear all
close all
clc
toolsPath = @(x) ...
    strcat('C:\Users\John\Documents\Astro\ASEN5050\tools\',x);
if ispc
    addpath(toolsPath(''))
end
% Cell array to track what functions are used, so they can be published
% later
global function_list;
function_list = {};
figWidth = 1120; % pixels
figHeight = 840; % pixels
CelestialConstants; % import useful constants

% Start with a regular Hohmann xfer for comparison. 
% Assume planets are in circular orbits. 
a_xfer = (Earth.a + Mars.a)/2;
v_earth_hoh = sqrt(2*Sun.mu/Earth.a - Sun.mu/a_xfer);
P_xfer = 2*pi*sqrt(a_xfer*a_xfer*a_xfer/Sun.mu);
X0 = [Earth.a;0;0;0;v_earth_hoh;0;0];

X0_polar = [Earth.a, 0, 0, v_earth_hoh/Earth.a];

if 0
pureTwoBodyPolar = @(t,X) polarProp(t,X,'No',0);
[t_array,X_array]=ode45(pureTwoBodyPolar,[0 P_xfer/2],X0_polar,options);
plot(X_array(:,1).*cos(X_array(:,2))/au2km,...
X_array(:,1).*sin(X_array(:,2))/au2km,'r')
end

planets = {Mars, Jupiter};
beta_array = [0.05 0.10 0.15];
line_styles = {'-','--',':'};
fig_handles = {};

for p_idx = 1:length(planets)
planet = planets{p_idx};

fprintf([planet.name ':\n'])

for idx = 1:length(beta_array)
fig_handles{p_idx,idx} = figure('Position', [0, 0, figWidth, figHeight]); % Transfer plot
plot(sind(1:360),cosd(1:360),'LineWidth',2); hold on; axis equal
plot(sind(1:360)*planet.a/au2km,cosd(1:360)*planet.a/au2km,'r','LineWidth',2);
plot(0,0,'ko','LineWidth',2);
end

fig_handles{p_idx,4} = figure('Position', [0, 0, figWidth, figHeight]); % Coning angle plot

% Start with a regular Hohmann xfer for comparison. 
% Assume planets are in circular orbits. 
a_xfer = (Earth.a + planet.a)/2;
v_earth = sqrt(Sun.mu/Earth.a);
v_earth_hoh = sqrt(2*Sun.mu/Earth.a - Sun.mu/a_xfer);
v_end_hoh = sqrt(2*Sun.mu/planet.a - Sun.mu/a_xfer);
v_planet = sqrt(Sun.mu/planet.a);
dv_hoh_tot = v_earth_hoh-v_earth+v_planet-v_end_hoh;
P_xfer = 2*pi*sqrt(a_xfer*a_xfer*a_xfer/Sun.mu);

fprintf(['Hohmann transfer for ' planet.name ':\n'])
fprintf('\tdV_final = %.3f km/s\n',v_end_hoh)
fprintf('\tdV_total = %.3f km/s\n',dv_hoh_tot)
fprintf('\ttransfer time = %.0f days\n',P_xfer/3600/24/2)
% X0 = [Earth.a;0;0;0;v_earth_hoh;0;0];

for b_idx = 1:length(beta_array)

beta = beta_array(b_idx);

% Set options on functions
pureTwoBodyPolar = @(t,X) polarProp(t,X,'Optimal',beta);
propagation_r_limit = @(t,X) detectDistLimit(t,X,planet.a);

% ODE45 options
tol=1e-12;
options=odeset('RelTol',tol,'AbsTol',[tol tol tol tol],...
'Events', propagation_r_limit);

X0 = [Earth.a, 0, 0, sqrt(Sun.mu/Earth.a/Earth.a/Earth.a)];
[t_array,X_array]=ode45(pureTwoBodyPolar,[0 P_xfer*6], X0,options);
figure(fig_handles{p_idx,b_idx})
plot(X_array(:,1).*cos(X_array(:,2))/au2km,...
X_array(:,1).*sin(X_array(:,2))/au2km,['k'],'LineWidth',2)

alpha_store = zeros(1,length(X_array));
% Xd_store = [];
for ii = 1:length(X_array)
    [Xd, alpha] = polarProp(0,X_array(ii,:)','Optimal',beta);
    alpha_store(ii) = alpha;
%     Xd_store(:,ii) = Xd;
    
end
figure(fig_handles{p_idx,4})
plot(X_array(:,2)/2/pi,alpha_store*180/pi, line_styles{b_idx},...
    'LineWidth',2); hold on;

v_mars = [0; sqrt(Sun.mu/Mars.a)]; % radial motion only
v_xfer_f = [X_array(end,3); (X_array(end,1)*X_array(end,4))];
dV_tot = abs(norm(v_mars-v_xfer_f));
transfer_time = t_array(end)/3600/24;

fprintf([planet.name ', beta = %.2f\n'],beta)
fprintf('\tdV_total = %.3f km/s\n',dV_tot)
fprintf('\ttransfer time = %.0f days\n',transfer_time)
fprintf('\tphase angle = %.0f degrees\n',X_array(end,2)*180/pi)
fprintf('\tintercept angle = %.0f degrees\n',atand(v_xfer_f(1)/v_xfer_f(2)))

end

for idx = 1:length(beta_array)
figure(fig_handles{p_idx,idx})
title(sprintf([planet.name ' Transfer, \\beta=%.2f'],beta_array(idx)))
xlabel('x (AU)')
ylabel('y (AU)')
legend('Earth orbit', [planet.name ' orbit'], 'Sun', 'Trajectory')
print(['Webpage/images/' planet.name num2str(idx)],'-dpng')
end

figure(fig_handles{p_idx,4})
title(['Sun Coning Angle for ' planet.name ' Transfer'])
xlabel('Orbit')
ylabel('\alpha (deg)')
legend(sprintf('\\beta=%.2f',beta_array(1)), ...
sprintf('\\beta=%.2f',beta_array(2)), ...
sprintf('\\beta=%.2f',beta_array(3)))
print(['Webpage/images/' planet.name 'ConeAngles'],'-dpng')
end

%% Jupiter
% 
% figure
% plot(sind(1:360),cosd(1:360)); hold on; axis equal
% plot(sind(1:360)*Jupiter.a/au2km,cosd(1:360)*Jupiter.a/au2km,'k');
% plot(0,0,'ko','LineWidth',2)