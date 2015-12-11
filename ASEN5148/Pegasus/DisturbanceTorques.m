%% Pegasus expected environmental torques
clear
close all
%% Constants
Re = 6378;
alt = 600;
r = Re+alt;
mu = 3.986e5;
A = 1.5; % Max area
d = 0.01; % Max CP/CG offset

%% GG
diff_I = 2;
T_gg = 3*mu/r/r/r*diff_I*pi/4

%% Radiation
P = 4.59e-6;
q = 0.3;
T_solar = P*A*d*(1+q)

%% drag
rho = 1.454e-13;
V = sqrt(mu/r)*1000;
Cd = 2.0;
D = 0.5*rho*V*V*Cd*A;
T_drag = D*d

%% magnetic
B0 = 3e-5;
L = pi/4; % 45 deg latitude (average)
% can simulate this over an orbit
B = B0*Re^3/r^3*sqrt(3*sin(L)^2+1);
M = 0.2; %0.2 from Brown
theta = pi/2;
T_mag = M*B*sin(theta)

%% s/c torques
% Actuators move the most at low beta angle
panel_area = 0.5; %m2
density = 3.5; %kg/m2
m = panel_area*density;
I = 1/12*m*panel_area; %only valid if panel is square...
P_over_2 = pi*sqrt(r^3/mu);
w_panels = -pi/P_over_2;
w_sc = pi/P_over_2;
% Multi-spin body equation from Schaub/Jenkins, no accel
T_sc = w_sc^2*(10+I+I) + 2*w_sc*w_panels*I




