%% ASEN 6080 Project Simulation -- John Clouse
% Estimation of a gravity-gradient microsat

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

[~, PV_state] = ode45(@two_body_state_dot, sim_tspan, [r_init; v_init],...
    PV_prop_opts.ode_opts, PV_prop_opts);

% Under gravity gradient + some disturbance torque

% A maneuver -- 180 deg slew about Z

m = 50; %kg
dim_length = 0.5; %m
dim_wd = sqrt(10)/10; %m
I11 = m/12*(dim_length^2+dim_wd^2);
I22 = m/12*(dim_length^2+dim_wd^2);
I33 = m/12*(dim_wd^2+dim_wd^2);