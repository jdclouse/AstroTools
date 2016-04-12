%% Jupiter mission by John Clouse
%% Initialize
CelestialConstants;
color_order = get(groot,'defaultAxesColorOrder');

% Set up the baseline events
% VEEJSN: Launch in 2015
JD_Launch = juliandate(2020,1,1); %July 1, 2015 00:00:00
JD_VGA = JD_Launch+100; % February 10, 1990 00:00:00
JD_VGA = juliandate([2020, 5, 30]);
JD_EGA1 = 2448235.5; % December 10, 1990 00:00:00
JD_EGA2 = 2448966.0; % December 9, 1992 12:00:00
JD_JOI = 2450164.0; % March 21, 1996 12:00:00

window_gran = 1;
num_resonant_years = 2;

% Anonymous function to print out the dates.
getDate = @(x_date) ...
    datestr(x_date-(floor(juliandate(date)) - datenum(date)),0);

hw_pub.figWidth = 1120; % pixels
hw_pub.figHeight = 840; % pixels
hw_pub.figPosn = [0, 0, hw_pub.figWidth, hw_pub.figHeight];
% Example: some_fig = figure('Position', hw_pub.figPosn);
hw_pub.lineWidth = 2; % pixels

%% Porkchop plots
% Launch to VGA
window = 0:window_gran:150; %was 150
Launch_dep = JD_Launch + window;
window = 0:window_gran:250;
VGA_arr = JD_VGA + window;
VGA_arr = juliandate([2020, 5, 30]):2:juliandate([2020, 10, 27]);

params1.fig_dim = hw_pub.figPosn;
params1.Sun = Sun;
params1.planet1 = Earth;
params1.planet2 = Venus;
params1.c3_countours = [10:1:13 15 18 25 40 50 60];
params1.v_inf_arr_countours = [1:1:9 11 17 29 45];
params1.v_inf_dep_countours = 11:1:23;
params1.TOF_countours = 50:50:500;
params1.day2sec = day2sec;
params1.show_c3 = true;
params1.show_v_inf_dep = false;
params1.show_v_inf_arr = true;
params1.show_tof = true;
params1.debug = false;

[ fh , output_Launch_VGA, legend_vec, legend_cells] = ...
    PorkchopPlot( Launch_dep, VGA_arr, params1);
figure(fh);
title('Launch to VGA');
%%
% VGA to EGA1
window = 0:window_gran:300;
EGA1_arr = JD_VGA + window + 150;%75 originally
params2.fig_dim = hw_pub.figPosn;
params2.Sun = Sun;
params2.planet1 = Venus;
params2.planet2 = Earth;
params2.c3_countours = 10:1:21;
params2.v_inf_arr_countours = [1:1:8 12 24];
params2.v_inf_dep_countours = [1:5 7 11 16 28];
params2.TOF_countours = 50:50:500;
params2.day2sec = day2sec;
params2.show_c3 = false;
params2.show_v_inf_dep = true;
params2.show_v_inf_arr = true;
params2.show_tof = true;
params2.debug = false;

[ fh , output_VGA_EGA1, legend_vec, legend_cells] = ...
    PorkchopPlot( VGA_arr, EGA1_arr, params2);
figure(fh);
xlabel('Departure Date')
title('VGA to EGA1');
%%
% EGA2 to JGA
EGA2_arr = EGA1_arr + 365*num_resonant_years;
% EGA2_arr = juliandate([2019, 12, 03]):2:juliandate([2020, 6, 01]);
window = 0:window_gran:(365*3);
window = 0:1:(365*3);
JGA_arr = EGA2_arr(1) + 365 + window;
% JGA_arr = juliandate([2020, 12, 03]):2:juliandate([2023, 4, 01]);
params3.fig_dim = hw_pub.figPosn;
params3.Sun = Sun;
params3.planet1 = Earth;
params3.planet2 = Jupiter;
params3.c3_countours = 10:1:21;
params3.v_inf_arr_countours = 1:1:10;
params3.v_inf_dep_countours = 10:2:18;
params3.TOF_countours = 300:200:2000;
params3.day2sec = day2sec;
params3.show_c3 = false;
params3.show_v_inf_dep = true;
params3.show_v_inf_arr = true;
params3.show_tof = true;
params3.debug = true;

[ fh , output_EGA2_JGA, legend_vec, legend_cells] = ...
    PorkchopPlot( EGA2_arr, JGA_arr, params3);
figure(fh);
xlabel('Departure Date')
title('EGA2 to JGA');

%% Patch together the trajectories using constraints
lambert_out = [output_Launch_VGA output_VGA_EGA1 output_EGA2_JGA];
% Initialize desired constraints
% can interp2 for more granularity
C3_max = 18; %km^2/s^2
launch_date_min = Launch_dep(1);%juliandate('9-Jan-2006');
launch_date_max = Launch_dep(end);%juliandate('10-Jan-2006');
V_inf_final_max = 6;
max_GA_diff = .1; %km/s
max_GA_diff = .15; %km/s

Earth_Venus.idx = 1;
Earth_Venus.route = 'long';
Earth_Venus.lambert = -1;
Venus_Earth.idx = 2;
Venus_Earth.route = 'long';
Venus_Earth.lambert = -1;

% Here we determine the constraints on the slice for the first segment
% Get the dates for valid C3
valid_c3_sw = (lambert_out(1).sw_c3_store < C3_max);
valid_c3_lw = (lambert_out(1).lw_c3_store < C3_max);
%apparently this is the element-wise logical AND 
valid_launch_lw = Launch_dep >= launch_date_min ...
    & Launch_dep < launch_date_max; 

get_val_c3 = @(leg) ['valid_c3_' leg.route(1) 'w'];


total_launch_constraints = ...
    eval(get_val_c3(Earth_Venus)) & repmat(valid_launch_lw,length(VGA_arr),1)';
launch_idx1 = find(sum(total_launch_constraints,2)>0, 1, 'first');
launch_idx2 = find(sum(total_launch_constraints,2)>0, 1, 'last');

fb_idx1 = 1;
fb_idx2 = length(EGA1_arr);

% Each segment forms a grid of information for each departure/arrival
% combination. For a gravity assist, this forms a 3D grid. However, the
% velocity difference must be small for pure-GA maneuvers. This section
% determines what the valid GAs are in this 3D grid. 
% for each launch date, see what has a valid GA
num_depart = length(Launch_dep);
num_VGA_window = length(VGA_arr);
num_EGA = length(EGA1_arr);
VGA_valid = zeros(num_depart, num_VGA_window, num_EGA);
VGA_vel_err_3d = nan(num_depart, num_VGA_window, num_EGA);

out_vel = @(leg) ['lambert_out(' num2str(leg.idx) ').' ...
    leg.route '_way_dv1_store'];
inc_vel = @(leg) ['lambert_out(' num2str(leg.idx) ').' ...
    leg.route '_way_dv2_store'];

incoming_v = eval(inc_vel(Earth_Venus));
outgoing_v = eval(out_vel(Venus_Earth));
for ii = launch_idx1:launch_idx2
    % 
    for jj = fb_idx1:fb_idx2
        GA_vel_err = abs(outgoing_v(:,jj) - incoming_v(ii,:)');
        valid_GA = ...
            GA_vel_err <= max_GA_diff;
        VGA_vel_err_3d(ii,:,jj) = GA_vel_err;
        VGA_valid(ii,:,jj) = valid_GA;
    end
end
%%
% 3:1 resonant orbit
% Need to keep the EGA windows the same. This will align the indices for
% EGA1 and EGA2.
% EGA2
% for simplicity, I had the windows the same. going to keep the indices the
% same as well, could offset half a day though...
num_JGA_window = length(JGA_arr);
ResoOrb_valid = zeros(num_VGA_window, num_EGA, num_JGA_window);
ResoOrb_vel_err = [];%nan(num_joi,1);
ResoOrb_vel_err_3d = nan(num_VGA_window, num_EGA, num_JGA_window);

% max_vinf_reso = 8;
% max_GA_diff_reso = .5;

Earth_Jupiter.idx = 3;
Earth_Jupiter.route = 'long';
Earth_Jupiter.lambert = -1;

incoming_v = eval(inc_vel(Venus_Earth));
outgoing_v = eval(out_vel(Earth_Jupiter));
for ii = 1:num_VGA_window
    % 
    for jj = 1:num_JGA_window;
        ResoOrb_vel_err = abs(outgoing_v(:,jj) - incoming_v(ii,:)');
        valid_RO = ...
            ResoOrb_vel_err <= max_GA_diff;
        ResoOrb_vel_err_3d(ii,:,jj) = ResoOrb_vel_err;
        ResoOrb_valid(ii,:,jj) = valid_RO;
    end
end
