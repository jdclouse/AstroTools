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
CelestialConstants; % import useful constants

W = 1361; %W/m2 = kg*m2/s3/m2
% W = 0; %W/m2 = kg*m2/s3/m2

% Spacecraft Specs from Wie
L_boom = 30;
h = L_boom*sin(pi/4)*2;
m_s = 40;
I = m_s*h*h/12;
% I = I + ms*2.5^2

m_p = 116; %kg
m = m_p+m_s;
J_s = I; %kg*m2
J_p = 20; %kg*m2
r = 0.88;
s = 0.94;
Bf = 0.79;
Bb = 0.55;
ef = 0.05;
eb=0.55;
rho_s = r*s;
rho_d = (Bf*r*(1-s)+(ef*Bf-eb*Bb)/(ef+eb))*3/2;
A_sail = h*h;
l = 2; %m
d = l*m_p/m;

SRP_Plot = figure;
dist_array = linspace(Earth.a,Mars.a)./au2km;
plot(dist_array, W/speed_of_light./(dist_array.*dist_array))
ylabel('Pressure (N/m)')
xlabel('Distance (AU)')

% The state vector will be X = [r; v; commanded_sun_angle];

% Anonymous function to calculate SRP
SRP = @(X) W/speed_of_light/(norm(X(1:3)/au2km));

% Calculate body acceleration
SRP_accel_body = @(X) SRP(X)*A_sail/m*[...
    (1+rho_s)*cos(X(7))*cos(X(7))+2/3*rho_d*cos(X(7));...
    (1-rho_s)*cos(X(7))*sin(X(7)); 0]; % return 3-vec
    
% Acceleration in inertial frame
SRP_accel_inrtl = @(X,t) ...
    Euler2DCM('3',getTrueAnom(X(1:3),X(4:6),Sun.mu,t)+X(7))...
    *SRP_accel_body(X); % return 3-vec

% Anonymous function to calculate 2-body accel
two_body = @(t,X) [X(4);X(5);X(6);...
    -Sun.mu*X(1)/norm(X(1:3))^3;...
    -Sun.mu*X(2)/norm(X(1:3))^3;...
    -Sun.mu*X(3)/norm(X(1:3))^3; 0]...
    + [0;0;0;SRP_accel_inrtl(X,t);0]; % return state-sized vec

% ODE45 options
tol=1e-12;
options=odeset('RelTol',tol,'AbsTol',[tol tol tol tol tol tol tol]);

% Start with a regular Hohmann xfer. 
% Assume planets are in circular orbits. 

a_xfer = (Earth.a + Mars.a)/2;
v_earth_hoh = sqrt(2*Sun.mu/Earth.a - Sun.mu/a_xfer);
P_xfer = 2*pi*sqrt(a_xfer*a_xfer*a_xfer/Sun.mu);
X0 = [Earth.a;0;0;0;v_earth_hoh;0;0];

fprintf('Pure Hohmann\n')
pureTwoBody = @(t,X) ProjectTwoBody(t,X,'No');
[t_array,X_array]=ode45(pureTwoBody,[0 P_xfer/2],X0,options);
fprintf('\n')

figure
plot(sind(1:360),cosd(1:360)); hold on; axis equal
plot(sind(1:360)*Mars.a/au2km,cosd(1:360)*Mars.a/au2km,'k');
plot(X_array(:,1)/au2km,X_array(:,2)/au2km,'r')

fprintf('Sun angle = 0\n')
twoBodySRP = @(t,X) ProjectTwoBody(t,X,'Yes');
X0 = [Earth.a;0;0;0;v_earth_hoh;0;0];
options=odeset('RelTol',tol,'AbsTol',[tol tol tol tol tol tol tol]);%,...
    %'Events', @detectMarsOrbit);
[t_array,X_array]=ode45(twoBodySRP,[0 P_xfer*2],X0,options);
fprintf('Done with ODE soln\n')
plot(X_array(:,1)/au2km,X_array(:,2)/au2km,'r')
fpa_store = [];
phase_store = [];
srp_accel_store = [];
% for ii = 1:length(X_array)
% %     r(ii) = norm(X_array(:,1:3))/au2km;
%     [~, srp_accel, fpa, phase_angle] = ProjectTwoBody(0,X_array(ii,:)','Optimal');
%     srp_accel_store = [srp_accel_store srp_accel(1:2)];
%     fpa_store(ii) = fpa;
%     phase_store(ii) = phase_angle;
% end
% figure
% plot(t_array,phase_store*180/pi)
% figure
% plot(t_array,fpa_store*180/pi)
% figure
% plot(t_array,srp_accel_store(1,:)); hold on
% plot(t_array,srp_accel_store(2,:)); 
%%
X0_polar = [Earth.a, 0, 0, v_earth_hoh/Earth.a];
pureTwoBodyPolar = @(t,X) polarProp(t,X,'No');
options=odeset('RelTol',tol,'AbsTol',[tol tol tol tol]);
[t_array,X_array]=ode45(pureTwoBodyPolar,[0 P_xfer/2],X0_polar,options);

figure
plot(sind(1:360),cosd(1:360)); hold on; axis equal
plot(sind(1:360)*Mars.a/au2km,cosd(1:360)*Mars.a/au2km,'k');
plot(X_array(:,1).*cos(X_array(:,2))/au2km,...
X_array(:,1).*sin(X_array(:,2))/au2km,'r')

pureTwoBodyPolar = @(t,X) polarProp(t,X,'Zero');
[t_array,X_array]=ode45(pureTwoBodyPolar,[0 P_xfer*2],X0_polar,options);
plot(X_array(:,1).*cos(X_array(:,2))/au2km,...
X_array(:,1).*sin(X_array(:,2))/au2km,'r')

% pureTwoBodyPolar = @(t,X) polarProp(t,X,'Optimal');
% [t_array,X_array]=ode45(pureTwoBodyPolar,[0 P_xfer*2],X0_polar,options);
% plot(X_array(:,1).*cos(X_array(:,2))/au2km,...
% X_array(:,1).*sin(X_array(:,2))/au2km,'r')

alpha_store = [];
for ii = 1:length(X_array)
%     r(ii) = norm(X_array(:,1:3))/au2km;
    [~, alpha] = polarProp(0,X_array(ii,:)','Optimal');
    alpha_store(ii) = alpha;
end
figure
plot(t_array,alpha_store*180/pi)
