%% Neptune mission by John Clouse
%% Initialize
CelestialConstants;
color_order = get(groot,'defaultAxesColorOrder');

% Set up the baseline events
% VEEJSN: Launch in 2015
JD_Launch = juliandate(2026,6,1); %July 1, 2015 00:00:00
JD_VGA = JD_Launch+50; % February 10, 1990 00:00:00
JD_EGA1 = 2448235.5; % December 10, 1990 00:00:00
JD_EGA2 = 2448966.0; % December 9, 1992 12:00:00
JD_JOI = 2450164.0; % March 21, 1996 12:00:00

window_gran = 1;
num_resonant_years = 2;

% Anonymous function to print out the dates.
getDate = @(x_date) ...
    datestr(x_date-(floor(juliandate(date)) - datenum(date)),0);

%% Porkchop plots
% Launch to VGA
window = 0:window_gran:150;
Launch_dep = JD_Launch + window;
window = 0:window_gran:250;
VGA_arr = JD_VGA + window;

params1.fig_dim = hw_pub.figPosn;
params1.Sun = Sun;
params1.planet1 = Earth;
params1.planet2 = Venus;
params1.c3_countours = [10:1:18 20:5:30 40 50 60];
params1.v_inf_arr_countours = [1:1:11 13 17 29 45];
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
EGA1_arr = VGA_arr + 75;
params2.fig_dim = hw_pub.figPosn;
params2.Sun = Sun;
params2.planet1 = Venus;
params2.planet2 = Earth;
params2.c3_countours = 10:1:21;
params2.v_inf_arr_countours = [2:2:12 17:4:21 29 37];
params2.v_inf_dep_countours = [1:11 12:4:28];
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
title('VGA to EGA1');
%%
% EGA2 to Mars
EGA2_arr = EGA1_arr + 365*num_resonant_years;
window = 0:window_gran:(365);
MGA_arr = EGA2_arr(1) + 160 + window;
params3.fig_dim = hw_pub.figPosn;
params3.Sun = Sun;
params3.planet1 = Earth;
params3.planet2 = Mars;
params3.c3_countours = 10:1:21;
params3.v_inf_arr_countours = 1:1:13;
params3.v_inf_dep_countours = 11:1:23;
params3.TOF_countours = 100:100:2000;
params3.day2sec = day2sec;
params3.show_c3 = false;
params3.show_v_inf_dep = true;
params3.show_v_inf_arr = true;
params3.show_tof = true;
params3.debug = false;

[ fh , output_EGA2_MGA, legend_vec, legend_cells] = ...
    PorkchopPlot( EGA2_arr, MGA_arr, params3);
figure(fh);
title('EGA2 to MGA');

%%
% MGA to JGA
window = 0:window_gran:(365*4);
JGA_arr = MGA_arr(1) + 365*3 + window;
params4.fig_dim = hw_pub.figPosn;
params4.Sun = Sun;
params4.planet1 = Mars;
params4.planet2 = Jupiter;
params4.c3_countours = 10:1:21;
params4.v_inf_arr_countours = 1:1:13;
params4.v_inf_dep_countours = 11:1:23;
params4.TOF_countours = 100:100:200;
params4.day2sec = day2sec;
params4.show_c3 = false;
params4.show_v_inf_dep = true;
params4.show_v_inf_arr = true;
params4.show_tof = true;
params4.debug = true;
params4.min_xfer_time = 300; % days
params4.lambert_tol = 1e-2;

fprintf('MGA to JGA Plot \n')
clock
tic
[ fh , output_MGA_JGA, legend_vec, legend_cells] = ...
    PorkchopPlot( MGA_arr, JGA_arr, params4);
% [ fh , output_JGA_SGA, legend_vec, legend_cells] = ...
%     PorkchopPlot( JGA_arr(69), SGA_arr(1047), params4);
toc
figure(fh);
title('MGA to JGA');

%% SGA to NOI
% NOIce.
window = 0:window_gran:(365*4);
UOI_arr = JGA_arr(1) + 365*6.5 + window;
% UOI_arr = juliandate(2027,12,23):UOI_arr(end);
params5.fig_dim = hw_pub.figPosn;
params5.Sun = Sun;
params5.planet1 = Jupiter;
params5.planet2 = Uranus;
params5.c3_countours = 10:1:21;
params5.v_inf_arr_countours = [1:1:15 16:2:30];
params5.v_inf_dep_countours = 3:2:25;
params5.TOF_countours = 100:1500:6000;
params5.day2sec = day2sec;
params5.show_c3 = false;
params5.show_v_inf_dep = true;
params5.show_v_inf_arr = true;
params5.show_tof = true;
params5.debug = true;
params5.min_xfer_time = 900; % days
params5.lambert_tol = 1e-2;
params5.ignore_longway = true;


fprintf('JGA to UOI Plot \n')
clock
tic
[ fh , output_JGA_UOI, legend_vec, legend_cells] = ...
    PorkchopPlot( JGA_arr, UOI_arr, params5);
% [ fh , output_JGA_SGA, legend_vec, legend_cells] = ...
%     PorkchopPlot( SGA_arr(1460), NOI_arr, params5);
toc
figure(fh);
title('JGA to UOI');

%% Patch together the trajectories using constraints
lambert_out = [output_Launch_VGA output_VGA_EGA1];% output_EGA2_JGA ...
%     output_JGA_SGA output_SGA_NOI];
% Initialize desired constraints
% can interp2 for more granularity
C3_max = 35; %km^2/s^2
launch_date_min = Launch_dep(1);%juliandate('9-Jan-2006');
launch_date_max = Launch_dep(end);%juliandate('10-Jan-2006');
V_inf_final_max = 14;
max_GA_diff = .4; %km/s

% Here we determine the constraints on the slice for the first segment
% Get the dates for valid C3
valid_c3_lw = (lambert_out(1).lw_c3_store < C3_max);
valid_c3_sw = (lambert_out(1).sw_c3_store < C3_max);
%apparently this is the element-wise logical AND 
valid_launch_lw = Launch_dep >= launch_date_min ...
    & Launch_dep < launch_date_max; 

total_launch_constraints = ...
    valid_c3_sw & repmat(valid_launch_lw,length(EGA1_arr),1)';
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

for ii = launch_idx1:launch_idx2
    % 
    for jj = fb_idx1:fb_idx2
        GA_vel_err = abs(lambert_out(2).short_way_dv1_store(:,jj) ...
            - lambert_out(1).long_way_dv2_store(ii,:)');
        valid_GA = ...
            GA_vel_err <= max_GA_diff;
        VGA_vel_err_3d(ii,:,jj) = GA_vel_err;
        VGA_valid(ii,:,jj) = valid_GA;
    end
end
%%

% 2:1 resonant orbit
% Need to keep the EGA windows the same. This will align the indices for
% EGA1 and EGA2.
% EGA2
% for simplicity, I had the windows the same. going to keep the indices the
% same as well, could offset half a day though...
num_MGA_window = length(MGA_arr);
ResoOrb_valid = zeros(num_depart, num_VGA_window, num_EGA);
ResoOrb_vel_err = [];%nan(num_joi,1);
ResoOrb_vel_err_3d = nan(num_VGA_window, num_EGA);

for ii = 1:num_VGA_window
    % 
    for jj = 1:num_MGA_window;
        ResoOrb_vel_err = ...
            abs(lambert_out(3).long_way_dv1_store(:,jj) ...
            - lambert_out(2).short_way_dv2_store(ii,:)');
        valid_RO = ...
            ResoOrb_vel_err <= max_GA_diff;
        ResoOrb_vel_err_3d(ii,:,jj) = ResoOrb_vel_err;
        ResoOrb_valid(ii,:,jj) = valid_RO;
    end
end
%%
% JGA
num_JGA_window = length(MGA_arr);
num_SGA = length(JGA_arr);
JGA_valid = int8(zeros(num_EGA, num_JGA_window, num_SGA));
JGA_vel_err_3d = nan(num_EGA, num_JGA_window, num_SGA);

for ii = 1:num_EGA
    % 
    for jj = 1:num_SGA
        GA_vel_err = abs(lambert_out(4).short_way_dv1_store(:,jj) ...
            - lambert_out(3).long_way_dv2_store(ii,:)');
        valid_GA = ...
            GA_vel_err <= max_GA_diff;
        JGA_vel_err_3d(ii,:,jj) = GA_vel_err;
        JGA_valid(ii,:,jj) = valid_GA;
    end
end


% SGA
num_NOI = length(UOI_arr);
SGA_valid = int8(zeros(num_JGA_window, num_SGA, num_NOI));
% SGA_vel_err_3d = nan(num_JGA_window, num_SGA, num_NOI);

for ii = 1:num_JGA_window
    % 
    for jj = 1:num_NOI
        GA_vel_err = abs(lambert_out(5).short_way_dv1_store(:,jj) ...
            - lambert_out(4).short_way_dv2_store(ii,:)');
        valid_GA = ...
            GA_vel_err <= max_GA_diff;
%         SGA_vel_err_3d(ii,:,jj) = GA_vel_err;
        SGA_valid(ii,:,jj) = valid_GA;
    end
end

%%
max_err = 1.5;
best_traj = zeros(6,1);
best_err = max_err;
lowest_vf_traj = zeros(6,1);
lowest_vf = V_inf_final_max;
lowest_vf_err = max_err;
lowest_C3_traj = zeros(6,1);
lowest_C3 = C3_max;
lowest_C3_err = max_err;
for launch_idx = launch_idx1:launch_idx2
    for VGA_idx = 1:num_VGA_window
        for EGA_idx = 1:num_EGA
            if VGA_valid(launch_idx, VGA_idx, EGA_idx) == 0
                continue
            end                
            for JGA_idx = 1:num_JGA_window
                if ResoOrb_valid(VGA_idx, EGA_idx, JGA_idx) == 0
                    continue
                end                
                for SGA_idx = 1:num_SGA
                    if JGA_valid(EGA_idx, JGA_idx, SGA_idx) == 0
                        continue
                    end                
                    for NOI_idx = 1:num_NOI
                        if SGA_valid(JGA_idx, SGA_idx, NOI_idx) == 0
                            continue
                        end  
                        % Congrats. you made it!
                        traj_error = ...
                            +VGA_vel_err_3d(launch_idx, VGA_idx, EGA_idx)...
                            +ResoOrb_vel_err_3d(VGA_idx, EGA_idx, JGA_idx)...
                            +JGA_vel_err_3d(EGA_idx, JGA_idx, SGA_idx)...
                            +abs(lambert_out(5).short_way_dv1_store(SGA_idx,NOI_idx) ...
                                - lambert_out(4).short_way_dv2_store(JGA_idx,SGA_idx)');
                        
                        if traj_error < max_err &&...
                                lambert_out(1).lw_c3_store(launch_idx,VGA_idx) < C3_max &&...
                                lambert_out(5).short_way_dv1_store(SGA_idx,NOI_idx) < V_inf_final_max
                            if traj_error < best_err
                                best_err = traj_error;
                                best_traj = [launch_idx;
                                    VGA_idx; EGA_idx; JGA_idx; SGA_idx; NOI_idx];
                            end
                            if lambert_out(1).lw_c3_store(launch_idx,VGA_idx) < lowest_C3_err
                                lowest_C3_err = lambert_out(1).lw_c3_store(launch_idx,VGA_idx);
                                lowest_C3_traj = [launch_idx;
                                    VGA_idx; EGA_idx; JGA_idx; SGA_idx; NOI_idx];
                            end
                            if lambert_out(5).short_way_dv1_store(SGA_idx,NOI_idx) < lowest_vf_err
                                lowest_vf_err = lambert_out(5).short_way_dv1_store(SGA_idx,NOI_idx);
                                lowest_vf_traj = [launch_idx;
                                    VGA_idx; EGA_idx; JGA_idx; SGA_idx; NOI_idx];
                            end
                        end
                    end
                end
            end
        end
    end
end
