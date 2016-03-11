%% Galileo Reconstruction by John Clouse
%% Initialize
CelestialConstants;

% Set up the baseline events
JD_Launch = 2447807.5; %October 8, 1989 00:00:00
JD_VGA = 2447932.5; % February 10, 1990 00:00:00
JD_EGA1 = 2448235.5; % December 10, 1990 00:00:00
JD_EGA2 = 2448966.0; % December 9, 1992 12:00:00
JD_JOI = 2450164.0; % March 21, 1996 12:00:00

%% Porkchop plots
% Launch to VGA
offset = 60;
window = -offset:0.5:offset;
Launch_dep = JD_Launch + window;
VGA_arr = JD_VGA + window;

params1.fig_dim = hw_pub.figPosn;
params1.Sun = Sun;
params1.planet1 = Earth;
params1.planet2 = Venus;
params1.c3_countours = 10:1:21;
params1.v_inf_arr_countours = 1:1:13;
params1.v_inf_dep_countours = 11:1:23;
params1.TOF_countours = 50:50:500;
params1.day2sec = day2sec;
params1.show_c3 = true;
params1.show_v_inf_dep = false;
params1.show_v_inf_arr = true;
params1.show_tof = true;
params1.debug = false;

[ fh , output_Launch_VGA] = PorkchopPlot( Launch_dep, VGA_arr, ...
    params1);
figure(fh);
title('Launch to VGA');

% VGA to EGA1
EGA1_window = JD_EGA1 + window;

params2 = params1;
params2.show_c3 = false;
params2.show_v_inf_dep = true;
params2.planet1 = Venus;
params2.planet2 = Earth;
params2.v_inf_dep_countours = [1:11 15 19 23];

[ fh , output_Launch_VGA] = PorkchopPlot( VGA_arr, EGA1_window, ...
    params2);
figure(fh);
title('VGA to EGA1');
