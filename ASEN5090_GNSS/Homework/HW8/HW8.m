%% HW 8
% John Clouse
% 
% read_GPSbroadcast, broadcast2xv, adjust year, and read_rinex functions 
% provided from class are used in this homework.

%% Initialize
clearvars -except function_list pub_opt 
close all
c = 2.99792458e8; %m/s
f_L1 = 1575.42; % MHz
f_L2 = 1227.60; % MHz

print_vec3 = @(x) ...
    fprintf(sprintf('%13.2f\n%13.2f\n%13.2f\n',x(1),x(2),x(3)));
print_vec4 = @(x) ...
    fprintf(sprintf('%13.2f\n%13.2f\n%13.2f\n%13.2f\n',x(1),x(2),x(3),x(4)));
%% Read the data files

eph = read_GPSbroadcast('brdc2550.14n');
obs_data = read_rinex_obs('test.14o');
P1 = 6;
P2 = 7; 
C1 = 8;
approx_rx_pos = [-1248596.2520 -4819428.2840  3976506.0340]'; % m
GPS_Week = eph(1,19);
GPS_TOD = [1 03 00];
TOW = eph(1,20)+GPS_TOD(1)*3600 + GPS_TOD(2)*60 + GPS_TOD(3);
fprintf('RINEX location:\n')
print_vec3(approx_rx_pos)
fprintf('\n')

%% Observations at the given epoch
obs_at_epoch = obs_data.data(obs_data.data(:,2) == TOW, :); 
prn_list = obs_at_epoch(:,3);
num_sats = length(prn_list);

% Ionosphere
iono_free_pseudorange = ...
    (f_L1*f_L1*obs_at_epoch(:,P1) - f_L2*f_L2*obs_at_epoch(:,P2))./ ...
    (f_L1*f_L1-f_L2*f_L2);

% Geometric range, satellite clock correction, relativity correction
geo_range = zeros(num_sats,1);
relativityCorr = zeros(num_sats,1);
satClkCorr = zeros(num_sats,1);
r_sat = zeros(num_sats,3);
for ii = 1:num_sats
    [~, geo_range(ii), tmp] = compute_range(eph, prn_list(ii), TOW, approx_rx_pos);
    r_sat(ii,:) = tmp';
    [~,~,~,relativityCorr(ii),satClkCorr(ii)] = ...
        broadcast2xv(eph,[GPS_Week TOW],prn_list(ii));
end

% Azimuth, Elevation
GPSvec = [2014 09 12 1 03 00];
% navfilename = generate_GPSyuma_name(GPSvec);
[navfilename,statusflag] = download_GPSyuma(GPSvec);
durationhrs = 1;
dt_sec = 3601;
ant_enu = [0 0 1];
mask_min = 0;  % deg
mask_max = 90; % deg
[latgd, lon, alt] = ECEF2ellipsoidal(approx_rx_pos);

[time_wntow, GPSdata] = ...
    ASEN5090_GPSvis(navfilename, 1, GPSvec,...
    durationhrs, dt_sec, latgd*180/pi, lon*180/pi, alt,...
    mask_min, mask_max, mask_min, ant_enu, 0, []);
% rearrange the data to be in the same order as prn list
prn_el = zeros(num_sats,1);
prn_az = zeros(num_sats,1);
for ii = 1:num_sats
    prn_el(ii) = GPSdata.topo_el(prn_list(ii));
    prn_az(ii) = GPSdata.topo_az(prn_list(ii));
end

T = table(prn_list, iono_free_pseudorange, geo_range, satClkCorr, ...
    relativityCorr, prn_el, prn_az);
fprintf('PRNs and corrections (distances in meters, angles in degrees):\n')
disp(T)
fprintf('\n\n\n\n\n')

%% Observation/Geometry matrix
A = ones(num_sats, 4); % x, y, z, rx time error
for ii = 1:num_sats
    A(ii,1:3) = -(r_sat(ii,:)-approx_rx_pos')./geo_range(ii);
end
A

%% Prefit residuals
% Plotting the residuals and Lecture 15's simple Tropo model vs. elevation
dy = iono_free_pseudorange - geo_range + satClkCorr - relativityCorr;

plot(prn_el, dy, '.')
hold on
T_el = 10:90;
plot(T_el, 2.5./sin(T_el*pi/180), 'r')
ylabel('(m)')
xlabel('Elevation (degrees)')
legend('Pre-fit residuals','Tropospheric delay')

%% Least squares solution for position deviation
% Find the position deviation from the RINEX location, dx

est_deviation = (A'*A)\A'*dy;
est_deviation_with_T = (A'*A)\A'*(dy - 2.5./sin(prn_el*pi/180));

fprintf('Estimated Deviation (m):\n')
print_vec4(est_deviation)
fprintf('Estimated Deviation with Troposphere model (m):\n')
print_vec4(est_deviation_with_T)
fprintf('\n')

%% Starting from incorrect location
% Choosing the center of the earth as the initial guess.
approx_rx_pos = [0 0 0]'; % m
b = 0;
for iter = 1:5
    for ii = 1:num_sats
        [~, geo_range(ii), tmp] = ...
            compute_range(eph, prn_list(ii), TOW, approx_rx_pos);
        r_sat(ii,:) = tmp';
        A(ii,1:3) = -(r_sat(ii,:)-approx_rx_pos')./geo_range(ii);
    end
    dy = iono_free_pseudorange - geo_range + satClkCorr - relativityCorr ;
    est_deviation = (A'*A)\A'*dy;
    approx_rx_pos = approx_rx_pos + est_deviation(1:3);
    b = est_deviation(4) + b;
    fprintf('Iteration %d\n',iter)
    fprintf('Estimated Deviation (m):\n')
    print_vec4(est_deviation)
    fprintf('New Position (m, ECF):\n')
    print_vec3(approx_rx_pos)
    fprintf('\n')
end