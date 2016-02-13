%% John Clouse IMD Lab3 
%% Initialize
clearvars -except hw_pub function_list

CelestialConstants

params.fig_dim = hw_pub.figPosn;
params.Sun = Sun;
params.planet1 = Earth;
params.planet2 = Mars;
params.c3_countours = [16 17 19 21];
params.v_inf_countours = [2.5 3 4 5 7.5];
params.TOF_countours = 50:50:500;
params.day2sec = day2sec;

% verification that the plot matches what's in the lab.
JD_dep1 = 2453525;
JD_arr1 = 2453705;

DepartureDates = (0:140) + JD_dep1;
ArrivalDates = (0:450) + JD_arr1;

PorkchopPlot( DepartureDates, ArrivalDates, params);
title('Verification example');

%% Problem 2
JD_dep1 = 2458200;
JD_arr1 = 2458350;

DepartureDates = JD_dep1:2458320;
ArrivalDates = JD_arr1:2458600;

params.c3_countours = [7 8 10 12 14 16 17 19 21];
params.v_inf_countours = [2.5 3 3.5 4 4.5 5 7.5];

PorkchopPlot( DepartureDates, ArrivalDates, params);
title('Problem 2');

%% Problem 3
JD_dep1 = 2457389;
JD_arr1 = 2457570;

DepartureDates = JD_dep1:2457509;
ArrivalDates = JD_arr1:2457790;

params.c3_countours = [7 8 10 12 14 16 17 19 21];
params.v_inf_countours = [2.5 3 3.5 4 4.5 5 7.5];

PorkchopPlot( DepartureDates, ArrivalDates, params);
title('Problem 3');
