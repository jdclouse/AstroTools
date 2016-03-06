clear

theta0 = 0;
theta_dot = 7.29211585275553e-005;
JD_init = 2456296.25;
mu_sun = 132712440017.987; % km3/s2
solar_flux = 1357; %W/m2 @ 1 AU
A_m_ratio =  0.01; % m2/kg
Re = 6378.1363;

SOI = 925000; %km

c = 299792.458; %k/s
mu_earth = 3.98600432896939e5; %km3/s2

EMO2EME_theta = 23.4393; %deg
EMO2EME = [1 0 0; 0 cosd(EMO2EME_theta) -sind(EMO2EME_theta); ...
    0 sind(EMO2EME_theta) cosd(EMO2EME_theta)];

au2km = 149597870.700; % km/AU

filename = 'Interplanetary_Obs.txt';
delimiterIn = ',';
headerlinesIn = 2;
y = importdata(filename,delimiterIn,headerlinesIn);
Obs=y.data;

ObsMassaged = zeros(length(Obs), 5);
ObsMassaged(:,2) = Obs(:,1);
for ii = 1:length(Obs)
    for jj = 1:3
        if ~isnan(Obs(ii,jj+1))
            ObsMassaged(ii,1) = jj;
            ObsMassaged(ii,3) = Obs(ii,jj+1);
            ObsMassaged(ii,4) = Obs(ii,jj+4);
            break;
        end
    end
end

% i. Canberra Station (DSS 34) at latitude = -35.398333?
% , longitude = 148.981944?
% and altitude = 0.691750 km
% ii. Madrid Station (DSS 65) at latitude = 40.427222?
% , longitude = 4.250556?
% and altitude = 0.834539 km
% 1
% iii. Goldstone Station (DSS 13) at latitude = 35.247164?
% , longitude = 243.205?
% and altitude = 1.07114904 km

site(1).name = 'Canberra Station (DSS 34)';
site(1).id = 1;
site(1).lat_lon_alt = [-35.398333; 148.981944;0.691750]; % geodetic, deg, km
site(1).r = lla2ecef(site(1).lat_lon_alt', 0, Re); % m
site(2).name = 'Madrid Station (DSS 65)';
site(2).id = 2;
site(2).lat_lon_alt = [40.427222; 4.250556; 0.834539]; % geodetic, deg, km
site(2).r = lla2ecef(site(2).lat_lon_alt', 0, Re); % m
site(3).name = 'Goldstone Station (DSS 13)';
site(3).id = 3;
site(3).lat_lon_alt = [35.247164; 243.205; 1.07114904]; % geodetic, deg, km
site(3).r = lla2ecef(site(3).lat_lon_alt', 0, Re); % m
num_sites = length(site);

filter_params;
propagator_opts.J2.use = 0;
propagator_opts.J3.use = 0;
propagator_opts.mu = mu_earth; %km3/s2
filter_opts.use_EKF = 0;
filter_opts.use_SNC = 0;

propagator_opts.Earth.Meeus.J200.L = [100.466449 35999.3728519 -0.00000568 0.0]; %deg
propagator_opts.Earth.Meeus.J200.a = 1.000001018*au2km; %km
propagator_opts.Earth.Meeus.J200.e = [0.01670862 -0.000042037 -0.0000001236 0.00000000004];
propagator_opts.Earth.Meeus.J200.i = [0 0.0130546 -0.00000931 -0.000000034]; % deg
propagator_opts.Earth.Meeus.J200.RAAN = [174.873174 -0.2410908 0.00004067 -0.000001327]; %deg
propagator_opts.Earth.Meeus.J200.Pi = [102.937348 0.3225557 0.00015026 0.000000478]; %deg
propagator_opts.epoch = JD_init;
propagator_opts.Sun.mu = mu_sun; % km3/s2
propagator_opts.solar_flux = solar_flux; %W/m2 @ 1 AU
propagator_opts.A_m_ratio =  A_m_ratio; % m2/kg
propagator_opts.c = c; %km/s
propagator_opts.au2km = au2km; %km
propagator_opts.EMO2EME = EMO2EME;

propagator_opts.OD.use = 1;
propagator_opts.OD.state_len = 7;
propagator_opts.OD.A_mat_handle = @A_state_rvCr;
filter_opts.H_tilda_handle = @H_tilda_state_rvCr;
propagator_opts.OD.A_params.mu = propagator_opts.mu;
propagator_opts.OD.A_params.Re = Re;
filter_opts.important_block = [7 7]; %rows, cols
propagator_opts.OD.A_params.important_block = filter_opts.important_block;
propagator_opts.OD.A_params.state_len = propagator_opts.OD.state_len;
STM_i = eye(propagator_opts.OD.state_len);
filter_opts.integrate_ref_state = 0;

% Filter Options
filter_opts.propagator_opts = propagator_opts;
filter_opts.ode_opts = odeset('RelTol', 1e-12, 'AbsTol', 1e-20);

filter_opts.theta_dot = theta_dot; %rad/s % TODO clean this!

state_ap = [-274096790.0 ; %km
            -92859240.0 ;
            -40199490.0 ;
            32.67 ; %km/s
            -8.94 ;
            -3.88 ;
            1.2];

P = diag([100 100 100 0.1 0.1 0.1 0.1]);
P = P.*P;

%% First iteration
% ref state
[~,X] = ode45(@flyby_two_body_state_dot, ObsMassaged(:,2), [state_ap; reshape(eye(7),49,1)], ...
        filter_opts.ode_opts, filter_opts.propagator_opts);
    
filter_opts.ref_state = X;
output_1 = SRIF(state_ap, P, ObsMassaged, filter_opts);
% figure
% plot3(X(:,1),X(:,2),X(:,3))
% xlabel('x (km)')
% ylabel('y (km)')
% zlabel('z (km)')
% axis equal

%% 2nd iteration
STM_accum = reshape(filter_opts.ref_state(end,7+1:end),...
    filter_opts.important_block(1), filter_opts.important_block(2));
x0_est = STM_accum\output_1.x_est_store(:,end);
iter2_state_ap = X(1,1:7)'+x0_est;
iter2_P = STM_accum\output_1.final_P/(STM_accum');

% ref state
[~,X] = ode45(@flyby_two_body_state_dot, ObsMassaged(:,2), [iter2_state_ap; reshape(eye(7),49,1)], ...
        filter_opts.ode_opts, filter_opts.propagator_opts);
  
filter_opts.ref_state = X;  
output_2 = SRIF(iter2_state_ap, P, ObsMassaged, filter_opts);

%% 3rd iteration
STM_accum = reshape(filter_opts.ref_state(end,7+1:end),...
    filter_opts.important_block(1), filter_opts.important_block(2));
x0_est = STM_accum\output_2.x_est_store(:,end);
iter3_state_ap = X(1,1:7)'+x0_est;
iter3_P = STM_accum\output_2.final_P/(STM_accum');

% ref state
[~,X] = ode45(@flyby_two_body_state_dot, ObsMassaged(:,2), [iter3_state_ap; reshape(eye(7),49,1)], ...
        filter_opts.ode_opts, filter_opts.propagator_opts);
  
filter_opts.ref_state = X;  
output_3 = SRIF(iter3_state_ap, P, ObsMassaged, filter_opts);

%% B-plane target

ode_opts = odeset('RelTol', 1e-12, 'AbsTol', 1e-20,...
    'Events',@stop_int);

[T, X_to_SOI] = ode45(@flyby_two_body_state_dot, ...
    [ObsMassaged(end,2), ObsMassaged(end,2)+100*86400], ...
    [output_3.state_store(:,end); reshape(eye(7),49,1)], ...
        ode_opts, filter_opts.propagator_opts);
    
% [~,~,~,~,~,~] = ...
%     cart2OE(X_to_SOI(end,1:3),...
%     X_to_SOI(end,4:6),mu_earth);
v_inf = norm(X_to_SOI(end,4:6));

f = acosh(1+ v_inf*v_inf/mu_earth* norm(X_to_SOI(end,1:3)));

LTOF = mu_earth/v_inf/v_inf/v_inf*(sinh(f)-f);

S_hat = X_to_SOI(end,4:6)/v_inf;
% T_hat
% The B Plane. note that without a V_out, this isn't very correct except
% the B Plane definition (so keep only that).
[~, ~, B_plane, ~, ~] = ...
    BPlaneTarget(X_to_SOI(end,4:6)', X_to_SOI(end,4:6)', mu_earth);

ode_opts = odeset('RelTol', 1e-12, 'AbsTol', 1e-20);
[~, X_to_BPlane] = ode45(@flyby_two_body_state_dot, ...
    [ObsMassaged(end,2), T(end)+LTOF], ...
    [output_3.state_store(:,end); reshape(eye(7),49,1)], ...
        ode_opts, filter_opts.propagator_opts);

% aim point: r_inf projected on the B plane
BT = dot(X_to_SOI(end,1:3)',B_plane(:,2));
BR = dot(X_to_SOI(end,1:3)',B_plane(:,3));

STM_DCO_intercept = reshape(X(end,8:end)',7,7);
P_intercept = STM_DCO_intercept*output_3.final_P*STM_DCO_intercept';
P_int_Bplane = B_plane*P_intercept(1:3,1:3)*B_plane';
[U,D] = eig(P_int_Bplane(2:3,2:3));
ellip_tilt = atan2(U(2,1),-U(1,1));

a=3*sqrt(max(max(D))); % horizontal radius
b=3*sqrt(min(max(D))); % vertical radius
x0=BT; % x0,y0 ellipse centre coordinates
y0=BR;
t=-pi:0.01:pi;
x=a*cos(t);
y=b*sin(t);
ell_rot = [cos(ellip_tilt) -sin(ellip_tilt); sin(ellip_tilt) cos(ellip_tilt)];
coords_prime = [];
for ii = 1:length(t)
    coords_prime(:,ii) = ell_rot*[x(ii); y(ii)];
end
figure
plot(x0+coords_prime(1,:),y0+coords_prime(2,:))
axis equal

%% plots
for ii = 1:7
figure
hold on
plot(output_1.state_store(ii,:))
plot(output_2.state_store(ii,:), 'r')
plot(output_3.state_store(ii,:), 'g')
end

fprintf('iter 1 x_est\n')
for ii = 1:7
    fprintf('%.5f\n',output_1.x_est_store(ii,end))
end
fprintf('\n')
fprintf('iter 1 x0_est\n')
for ii = 1:7
    fprintf('%.5f\n',x0_est(ii))
end
fprintf('\n')
fprintf('iter 3 X_final\n')
for ii = 1:7
    fprintf('%.5f\n',output_3.state_store(ii,end))
end