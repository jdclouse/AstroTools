%% John Clouse IMD HW 2, Problem 1
% 
%% Initialize
clearvars -except hw_pub function_list

CelestialConstants
JD = 2458239; %May 1, 2018
launch_date = '01-May-2018';
% JD = 2453587; %Aug 5, 2005
TOF_vec = 150:1:300; % Time of flight vector, days
num_elems = length(TOF_vec);
C3_limit = 50;

[re, ve] = MeeusEphemeris(Earth, JD, Sun);

type1_dv1_store = zeros(1,num_elems);
type1_dv2_store = zeros(1,num_elems);
type2_dv1_store = zeros(1,num_elems);
type2_dv2_store = zeros(1,num_elems);
cnt = 0;

%% Loop through all the days past the launch date
for TOF = TOF_vec
    cnt = cnt + 1;
    [rm, vm] = MeeusEphemeris(Mars, JD+TOF, Sun); %Add days to JD
    % Type 1
    [v1, v2] = lambert(re, rm, TOF*day2sec, 1, Sun);
    type1_dv1_store(cnt) = norm(v1-ve);
    type1_dv2_store(cnt) = norm(v2-vm);
    % Type 2
    [v1, v2] = lambert(re, rm, TOF*day2sec, -1, Sun);
    type2_dv1_store(cnt) = norm(v1-ve);
    type2_dv2_store(cnt) = norm(v2-vm);
end


%% Plot
figure('Position', hw_pub.figPosn)
subplot(2,1,1)
type1_C3 = type1_dv1_store.*type1_dv1_store;
plot(TOF_vec(type1_C3<C3_limit), type1_C3(type1_C3<C3_limit))
title('Type 1 trajectories')
ylabel('C3 (km^2/s^2)')
subplot(2,1,2)
plot(TOF_vec(type1_C3<C3_limit), type1_dv2_store(type1_C3<C3_limit))
ylabel('Arrival V_{inf} (km/s)')
xlabel(['Days past ' launch_date])

figure('Position', hw_pub.figPosn)
subplot(2,1,1)
type2_C3 = type2_dv1_store.*type2_dv1_store;
plot(TOF_vec(type2_C3<C3_limit), type2_C3(type2_C3<C3_limit))
title('Type 2 trajectories')
ylabel('C3 (km^2/s^2)')
subplot(2,1,2)
plot(TOF_vec(type2_C3<C3_limit), type2_dv2_store(type2_C3<C3_limit))
ylabel('Arrival V_{inf} (km/s)')
xlabel(['Days past ' launch_date])

%% Results
fprintf(['Minimum C3 for Type I xfer: %.1f km^2/s^2, arrival ' ...
    datestr(datenum(launch_date)+TOF_vec(type1_C3 == min(type1_C3))) '\n'], ...
    min(type1_C3))
fprintf(['Minimum V_inf for Type I xfer: %.3f km/s, arrival ' ...
    datestr(datenum(launch_date)+TOF_vec(type1_dv2_store == min(type1_dv2_store))) '\n'], ...
    min(type1_dv2_store))
fprintf(['Minimum C3 for Type II xfer: %.1f km^2/s^2, arrival ' ...
    datestr(datenum(launch_date)+TOF_vec(type2_C3 == min(type2_C3))) '\n'], ...
    min(type2_C3))
fprintf(['Minimum V_inf for Type II xfer: %.3f km/s, arrival ' ...
    datestr(datenum(launch_date)+TOF_vec(type2_dv2_store == min(type2_dv2_store))) '\n'], ...
    min(type2_dv2_store))