%% ASEN 6080 Project Simulation -- John Clouse
% Estimation of a gravity-gradient microsat
addpath('C:\Users\John\Documents\Astro\ASEN5050\tools');
addpath('C:\Users\John\Documents\Astro\ASEN6008_IMD\tools');
addpath('C:\Users\John\Documents\Astro\ASEN5010\AIAA_Software\Matlab Toolbox\UNIX-OSX');
addpath('C:\Users\John\Documents\Astro\ASEN5010\tools');
addpath('C:\Users\John\Documents\Astro\ASEN6080_SOD2\tools');

CelestialConstants;
quat_reset_dt = 1; % How often the quaternion is re-normalized.

% Circular orbit, point mass
orbit.a = Earth.R + 300; %km
orbit.e = 0;
orbit.i = 23*pi/180; %Should I make this sun sync?
orbit.w = 0;
orbit.RAAN = 0; 
orbit.f = 0;
orbit.P = 2*pi*sqrt(orbit.a^3/Earth.mu);
orbit.n = sqrt(Earth.mu/orbit.a^3);

sim_tspan = 0:quat_reset_dt:orbit.P;
[r_init, v_init] = OE2cart(orbit.a, orbit.e, orbit.i, orbit.RAAN, ...
    orbit.w, orbit.f, Earth.mu);

PV_prop_opts.mu = Earth.mu;
PV_prop_opts.J2.use = false;
PV_prop_opts.J3.use = false;
PV_prop_opts.drag.use = false;
PV_prop_opts.OD.use = false;
PV_prop_opts.ode_opts = odeset('RelTol', 1e-12, 'AbsTol', 1e-20);

% [~, PV_state] = ode45(@two_body_state_dot, sim_tspan, [r_init; v_init],...
%     PV_prop_opts.ode_opts, PV_prop_opts);

% Under gravity gradient + some disturbance torque

% A maneuver -- 180 deg slew about Z
% Spacecraft parameters
m = 50; %kg
dim_length = 0.5; %m
dim_wd = sqrt(10)/10; %m
I11 = m/12*(dim_length^2+dim_wd^2);
I22 = m/12*(dim_length^2+dim_wd^2);
I33 = m/12*(dim_wd^2+dim_wd^2);

% Initial attitude from LVLH
% att_prop_opts.PV_state = PV_state;
att_prop_opts.gravity_gradient.use = true;
att_prop_opts.mu = PV_prop_opts.mu;
att_prop_opts.I = [I11; I22; I33];

euler_angs = [0;5;-3]*pi/180;
% euler_angs = [0;0;-3]*pi/180;
% euler_angs = [0;0;0]*pi/180;
euler_seq = '321';
Qi_lvlh2body = Euler2EP(euler_seq, euler_angs);
Qi_inrtl2lvlh = C2EP(inrtl2lvlh(r_init, v_init));
Qi = addEP(Qi_inrtl2lvlh, Qi_lvlh2body);
ratei = [0; -orbit.n; 0.00655*pi/180];
ratei = [0; -orbit.n; 0];

combined_opts.PV_prop_opts = PV_prop_opts;
combined_opts.att_prop_opts = att_prop_opts;
sim_out = zeros(length(sim_tspan),13);
sim_out(1,:) = [r_init' v_init' Qi' ratei'];

% The simulation
for ii = 2:length(sim_tspan)
%     ii
    %sim_out(ii-1:ii,:)
    [~, iter_out] = ode45(@combined_state_dot, ...
        sim_tspan(ii-1:ii), sim_out(ii-1,:)',...
        PV_prop_opts.ode_opts, combined_opts);
    
    sim_out(ii,:) = iter_out(end,:);
    % renormalize the quaternion.
    sim_out(ii,7:10) = sim_out(ii,7:10)/norm(sim_out(ii,7:10));
end
% tic
%     [~, iter_out] = ode45(@combined_state_dot, ...
%         sim_tspan(1:10), sim_out(1,:)',...
%         PV_prop_opts.ode_opts, combined_opts);
%     toc

% Compute the euler angles for visualization
euler_angs = zeros(length(sim_tspan),3);
for ii = 1:length(sim_tspan)
    DCM = inrtl2lvlh(sim_out(ii,1:3)', sim_out(ii,4:6)');
    Q_inrtl2lvlh = C2EP(DCM);
    Q_lvlh2body = subEP(sim_out(ii,7:10), Q_inrtl2lvlh);
%     Q_lvlh2body = subEP(Q_inrtl2lvlh, sim_out(ii,7:10));
    euler_angs(ii,:) = EP2Euler321(Q_lvlh2body)';
end

figure; 
for ii = 1:3
subplot(3,1,ii)
plot(euler_angs(:,ii)*180/pi)
end

rate_meas_data = sim_out(:,11:13) + normrnd(0,1e-6,length(sim_tspan),3);
rate_meas_data = [sim_tspan' rate_meas_data];
