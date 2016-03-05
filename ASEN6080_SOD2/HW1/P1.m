%% HW1 Problem 1
% John Clouse
% Create orbit data set
%% Initialize
clear
global function_list;
function_list = {};
close all

mu = 3.986004418e5; %km3/s2
Re = 6378.137; %km
theta_dot = 7.2921158553e-5; %rad/s

theta = 0;

site(1).name = 'Boulder, CO';
site(1).id = 1;
site(1).lat_lon_alt = [40*pi/180;-105*pi/180;1.655]; % geodetic, rad, km
site(1).r = ellipsoidal2cart(site(1).lat_lon_alt(1),...
    site(1).lat_lon_alt(2),site(1).lat_lon_alt(3)) * 1e3; % m
site(2).name = 'Perth, Australia';
site(2).id = 2;
site(2).lat_lon_alt = [-32*pi/180;116*pi/180;.020]; % geodetic, rad, km
site(2).r = ellipsoidal2cart(site(2).lat_lon_alt(1),...
    site(2).lat_lon_alt(2),site(2).lat_lon_alt(3)) * 1e3; % m
site(3).name = 'Diego Garcia';
site(3).id = 3;
site(3).lat_lon_alt = [-7*pi/180;72.5*pi/180;.010]; % geodetic, rad, km
site(3).r = ellipsoidal2cart(site(3).lat_lon_alt(1),...
    site(3).lat_lon_alt(2),site(3).lat_lon_alt(3)) * 1e3; % m
site(4).name = 'KSC';
site(4).id = 4;
site(4).lat_lon_alt = [28.5*pi/180;-80.5*pi/180;.010]; % geodetic, rad, km
site(4).r = ellipsoidal2cart(site(4).lat_lon_alt(1),...
    site(4).lat_lon_alt(2),site(4).lat_lon_alt(3)) * 1e3; % m
site(5).name = 'Singapore';
site(5).id = 5;
site(5).lat_lon_alt = [(1+17/60)*pi/180;(103+50/60)*pi/180;.010]; % geodetic, rad, km
site(5).r = ellipsoidal2cart(site(5).lat_lon_alt(1),...
    site(5).lat_lon_alt(2),site(5).lat_lon_alt(3)) * 1e3; % m
site(6).name = 'Johannesburg, South Africa';
site(6).id = 6;
site(6).lat_lon_alt = [-26*pi/180;28*pi/180;.010]; % geodetic, rad, km
site(6).r = ellipsoidal2cart(site(6).lat_lon_alt(1),...
    site(6).lat_lon_alt(2),site(6).lat_lon_alt(3)) * 1e3; % m
num_sites = 6;

propagator_opts.OD.use = 0;
propagator_opts.drag.use = 0;
propagator_opts.J2.use = 1;
propagator_opts.J2.params.J2 = 0.0010826267;
propagator_opts.J2.params.mu = mu; 
propagator_opts.J2.params.Re = Re;
propagator_opts.J3.use = 1;
propagator_opts.J3.params.J3 = -2.5327e-6;
propagator_opts.J3.params.mu = mu; 
propagator_opts.J3.params.Re = Re;
propagator_opts.mu = mu; 

% Orbit characteristics
a = 10000; %km
e = 0.05;
i = 25*pi/180; %rad
RAAN = 210*pi/180; %rad
w = 35*pi/180; %rad
f = 0;
orb_period = 2*pi*sqrt(a*a*a/3.986e5);
% state_i = oe2cart([a,e,i,RAAN,w,f]);
state_i = zeros(6,1);
[state_i(1:3),state_i(4:6)] = OE2cart(a,e,i,RAAN,w,f,mu);

% Prop time
prop_time = 3600*24; % sec
meas_dt = 10; %sec
times = 0:meas_dt:prop_time;
num_pts = length(times);

ode_opts = odeset('RelTol', 1e-12, 'AbsTol', 1e-20);

% Propagate
[T,X] = ode45(@two_body_state_dot, times, state_i, ode_opts, propagator_opts);
J3_accel_store = zeros(3,num_pts);
%%
for ii = 1:num_pts
%     J3_accel_store(:,ii) = J2_accel( X(ii,1:3)', propagator_opts.J2.params );
    J3_accel_store(:,ii) = J3_accel( X(ii,1:3)', propagator_opts.J3.params );
end
figure
plot(times, J3_accel_store)

meas_store = []; % site, time, measurement
r_store = [];
el_mask = 10*pi/180; %rad
true_state = [];
for t_i = 1:num_pts
    for site_num = 1:6
        r_rel_ECEF = Euler2DCM('3',theta_dot*T(t_i))*X(t_i,1:3)' - site(site_num).r;
        r_rel_ENU = R_ECEF2ENU(site(site_num).lat_lon_alt(1),...
            site(site_num).lat_lon_alt(2))*r_rel_ECEF;
        if asin(r_rel_ENU(3)/norm(r_rel_ENU)) > el_mask
            meas_store = [meas_store; [site_num, T(t_i), ...
                1e3*compute_range_ECFsite(X(t_i,1:3)',site(site_num).r, ...
                theta_dot*T(t_i)),...
                1e3*compute_range_rate_ECFsite(X(t_i,1:6),site(site_num).r,... % m/s
                theta_dot*T(t_i), theta_dot), t_i]];
%                 1e3*norm(r_rel_ENU), ... % m
%                 compute_range_ECFsite(X(t_i,1:3)',site(site_num).r, ...
%                 theta_dot*T(t_i))]]; 
            true_state = [true_state, X(t_i,1:6)'];

        end
    end
    r_store(t_i) = norm(X(t_i,1:3));
end
% meas_store(:,3) = meas_store(:,3) + normrnd(0,10);
clean_meas_store(:,1) = meas_store(:,3);
clean_meas_store(:,2) = meas_store(:,4);
meas_store(:,3) = clean_meas_store(:,1) + normrnd(0,0.01,length(meas_store),1);
meas_store(:,4) = clean_meas_store(:,2) + normrnd(0,0.001,length(meas_store),1);

% Find 
figure
plot3(X(:,1),X(:,2),X(:,3))
hold on
plot3(cos(0:0.1:2*pi)*6378, zeros(length(0:0.1:2*pi)),sin(0:0.1:2*pi)*6378,'r')
plot3(zeros(length(0:0.1:2*pi)), cos(0:0.1:2*pi)*6378,sin(0:0.1:2*pi)*6378,'r')
plot3(state_i(1),state_i(2),state_i(3),'ro')
axis equal
% 
plot_arc_site = 2;
plot3(X(meas_store(meas_store(:,1)==plot_arc_site,5),1),...
    X(meas_store(meas_store(:,1)==plot_arc_site,5),2),...
    X(meas_store(meas_store(:,1)==plot_arc_site,5),3),'r*')
xlabel('X')
ylabel('Y')
zlabel('Z')

hold on
plot_arc_site = 1;
all_idx = 1:length(meas_store);
error = [];
for pass_idx = all_idx(meas_store(:,1)==plot_arc_site);
derp = X(meas_store(pass_idx,5),1:3)'/norm(X(meas_store(pass_idx,5),1:3))...
    *meas_store(pass_idx,3)*1e-3;
pos_inrtl = X(meas_store(pass_idx,5),1:3)';
site_inrtl = Euler2DCM('3',-theta_dot*meas_store(pass_idx,2))*site(meas_store(pass_idx,1)).r;
rel_vec = pos_inrtl - site_inrtl;
rel_vec = rel_vec/norm(rel_vec)*meas_store(pass_idx,3)*1e-3;
plot3([site_inrtl(1)],...
    [site_inrtl(2)],...
    [site_inrtl(3)],...
    'ko')
plot3([site_inrtl(1), site_inrtl(1)+rel_vec(1)],...
    [site_inrtl(2), site_inrtl(2)+rel_vec(2)],...
    [site_inrtl(3), site_inrtl(3)+rel_vec(3)],...
    'k')
error = [error; meas_store(pass_idx,3)*1e-3 - norm(pos_inrtl - site_inrtl)];
end
figure
plot(error)


figure
% plot(r_store)
subplot(6,1,1)
for selected_site = 1:6;
subplot(6,1,selected_site)
plot(meas_store(meas_store(:,1)==selected_site,2)/3600, ...
    meas_store(meas_store(:,1)==selected_site,3))
title([site(selected_site).name ' Range'])
ylabel('m')
end
xlabel('hr')

figure
% plot(r_store)
subplot(6,1,1)
for selected_site = 1:6;
subplot(6,1,selected_site)
plot(meas_store(meas_store(:,1)==selected_site,2)/3600, ...
    meas_store(meas_store(:,1)==selected_site,4))
title([site(selected_site).name ' Range Rate'])
ylabel('m/s')
end
xlabel('hr')