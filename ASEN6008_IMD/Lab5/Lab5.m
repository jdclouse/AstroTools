%% John Clouse IMD Lab 5 
% 
%% Initialize
clearvars -except hw_pub function_list

CelestialConstants

JD_launch = 2453755.29167;
JD_Jupiter = 2454159.73681; %JGA
JD_Pluto = 2457217.99931;

%% Planet P and V, Lambert solutions
[r_earth, v_earth] = MeeusEphemeris(Earth, JD_launch, Sun);
[r_jupiter, v_jupiter] = MeeusEphemeris(Jupiter, JD_Jupiter, Sun);
[r_pluto, v_pluto] = MeeusEphemeris(Pluto, JD_Pluto, Sun);

[v_earth_dep, v_jupiter_arr] = lambert(r_earth, r_jupiter, (JD_Jupiter - JD_launch)*day2sec, 1, Sun);
[v_jupiter_dep, v_pluto_arr] = lambert(r_jupiter, r_pluto, (JD_Pluto - JD_Jupiter)*day2sec, 1, Sun);

% Earth launch params
v_inf_earth_out = v_earth_dep - v_earth;
C3 = norm(v_inf_earth_out)^2;
RLA = atan2(v_inf_earth_out(2),v_inf_earth_out(1));
DLA = asin(v_inf_earth_out(3)/norm(v_inf_earth_out));

% JGA conditions
v_inf_JGA_in = v_jupiter_arr - v_jupiter;
v_inf_JGA_out = v_jupiter_dep - v_jupiter;

% Pluto flyby
v_inf_pluto_arr = v_pluto_arr - v_pluto;

%% Jupiter gravity assist
[b, B_hat, B_plane, turn_ang, r_JGA] = ...
    BPlaneTarget(v_inf_JGA_in, v_inf_JGA_out, Jupiter.mu);

R_jupiter = 71492; %km
h_JGA = r_JGA - R_jupiter;

helio_dv = norm(v_inf_JGA_out - v_inf_JGA_in)
helio_dv = norm(Euler2DCM('3',-turn_ang)*v_inf_JGA_in - v_inf_JGA_in)

%% Porkchop plots
params1.fig_dim = hw_pub.figPosn;
params1.Sun = Sun;
params1.planet1 = Earth;
params1.planet2 = Jupiter;
params1.c3_countours = [16 17 19 21];
params1.v_inf_arr_countours = 11:0.5:23;
params1.v_inf_dep_countours = 11:1:23;
params1.TOF_countours = 50:50:500;
params1.day2sec = day2sec;
params1.show_c3 = false;
params1.show_v_inf_dep = true;
params1.show_v_inf_arr = true;
params1.show_tof = true;

Earth_dep_dates = 2453714.5:1:2453794.5;
Jupiter_arr_dates = 2454129.5:1:2454239.5;
Pluto_arr_dates = 2456917.5:1:2457517.5;

[plot1, out1] = PorkchopPlot( Earth_dep_dates, Jupiter_arr_dates, ...
    params1);

params2 = params1;
params2.planet1 = Jupiter;
params2.planet2 = Pluto;
params2.TOF_countours = 50:50:4000;
params2.v_inf_arr_countours = 11:1:23;
params2.v_inf_dep_countours = 11:0.2:23;
[plot2, out2] = PorkchopPlot( Jupiter_arr_dates, Pluto_arr_dates, ...
    params2);

%% 
% Initialize desired constraints
% can interp2 for more granularity
% These are only short-ways... since it's for a long interplanetary mission
C3_max = 180; %km^2/s^2
launch_date_min = juliandate('9-Jan-2006');
launch_date_max = juliandate('10-Jan-2006');
V_flyby_final_max = 14.5;
max_GA_diff = .4; %km/s

% Here we determine the constraints on the slice for the first segment
% Get the dates for valid C3
valid_c3_sw = (out1.sw_c3_store < C3_max);
%apparently this is the element-wise logical AND 
valid_launch_sw = Earth_dep_dates >= launch_date_min ...
    & Earth_dep_dates < launch_date_max; 

total_launch_constraints = ...
    valid_c3_sw & repmat(valid_launch_sw,length(Jupiter_arr_dates),1)';
launch_idx1 = find(sum(total_launch_constraints,2)>0, 2, 'first');
launch_idx2 = find(sum(total_launch_constraints,2)>0, 2, 'last');

% Y of this plot is X of the next.
GA_idx1 = find(sum(total_launch_constraints,1)>0, 1, 'first');
GA_idx2 = find(sum(total_launch_constraints,1)>0, 1, 'last');

% Next determine the constraints on the slice for the last segment
valid_final_pass = out2.short_way_dv2_store < V_flyby_final_max;
fb_idx1 = find(sum(valid_final_pass,1)>0, 1, 'first');
fb_idx2 = find(sum(valid_final_pass,1)>0, 1, 'last');

% new_Earth_dates = Earth_dep_dates(launch_idx1):.1:Earth_dep_dates(launch_idx2+1);
% new_Jupiter_arr_dates = Jupiter_arr_dates(GA_idx1):.1:Jupiter_arr_dates(GA_idx2);
% new_Pluto_arr_dates = Pluto_arr_dates(fb_idx1):1:Pluto_arr_dates(fb_idx2);

% [plot3, out1_refined] = PorkchopPlot( new_Earth_dates, new_Jupiter_arr_dates, ...
%     params1);
% [plot4, out2_refined] = PorkchopPlot( new_Jupiter_arr_dates, new_Pluto_arr_dates, ...
%     params);

% Each segment forms a grid of information for each departure/arrival
% combination. For a gravity assist, this forms a 3D grid. However, the
% velocity difference must be small for pure-GA maneuvers. This section
% determines what the valid GAs are in this 3D grid. 
% for each launch date, see what has a valid GA
GA_valid = zeros(length(Earth_dep_dates),...
    111, length(Pluto_arr_dates));
GA_vel_err_3d = nan(length(Earth_dep_dates),...
    111, length(Pluto_arr_dates));
for ii = launch_idx1:launch_idx2
    % 
    for jj = fb_idx1:fb_idx2
        GA_vel_err = abs(out2.short_way_dv1_store(:,jj) ...
            - out1.short_way_dv2_store(ii,:)');
        valid_GA = ...
            GA_vel_err <= max_GA_diff;
        GA_vel_err_3d(ii,:,jj) = GA_vel_err;
        GA_valid(ii,:,jj) = valid_GA;
    end
end

% TODO: bail out if there isn't something close. either you have a bad
% range or not enough granularity in the dates.

% Those were velocities that matched. Now apply the valid dates to this.
GA_valid = GA_valid & repmat(total_launch_constraints, [1 1 length(Pluto_arr_dates)]);
% And the final passes.
GA_valid = GA_valid & repmat(reshape(valid_final_pass, 1, 111, 601), [length(Earth_dep_dates) 1 1]);

% Create 3D grids of the interesting parameters.
C3_X = GA_valid.*repmat(out1.sw_c3_store, [1 1 length(Pluto_arr_dates)]);
C3_X(C3_X == 0) = NaN;
viable_Vinf_pluto = GA_valid.*repmat(reshape(out2.short_way_dv2_store, 1, 111, 601), [length(Earth_dep_dates) 1 1]);
viable_Vinf_pluto(viable_Vinf_pluto == 0) = NaN;

% Lowest C3
[d3, PFB_date_idx] = min(C3_X, [], 3);
[d2, JGA_date_idx] = min(d3, [], 2);
[~, Launch_date_idx] = min(d2, [], 1);
Launch_date_idx
JGA_date_idx=JGA_date_idx(Launch_date_idx)
PFB_date_idx=PFB_date_idx(Launch_date_idx,JGA_date_idx)
C3_X(Launch_date_idx, JGA_date_idx, PFB_date_idx)

% Earliest to Pluto
for ii = 1:length(Pluto_arr_dates)
    if sum(sum(squeeze(GA_valid(:,:,ii))))
        break
    end
end
Earliest_Pluto = Pluto_arr_dates(ii);
jd_dn_delta = floor(juliandate(date)) - datenum(date);
datestr(Earliest_Pluto-jd_dn_delta,1)
getDate = @(x_date) datestr(x_date-(floor(juliandate(date)) - datenum(date)),1);
getDate(Earliest_Pluto)


% Least velocity
[d3, PFB_date_idx] = min(viable_Vinf_pluto, [], 3);
[d2, JGA_date_idx] = min(d3, [], 2);
[~, Launch_date_idx] = min(d2, [], 1);
Launch_date_idx
JGA_date_idx=JGA_date_idx(Launch_date_idx)
PFB_date_idx=PFB_date_idx(Launch_date_idx,JGA_date_idx)
viable_Vinf_pluto(Launch_date_idx, JGA_date_idx, PFB_date_idx)

% best trajectory
% I'll do it to minimize the error in the V_infinities
[d3, PFB_date_idx] = min(GA_vel_err_3d, [], 3);
[d2, JGA_date_idx] = min(d3, [], 2);
[~, Launch_date_idx] = min(d2, [], 1);
Launch_date_idx
JGA_date_idx=JGA_date_idx(Launch_date_idx)
PFB_date_idx=PFB_date_idx(Launch_date_idx,JGA_date_idx)

getDate(Pluto_arr_dates(PFB_date_idx))
getDate(Jupiter_arr_dates(JGA_date_idx))
getDate(Earth_dep_dates(Launch_date_idx))
C3_X(Launch_date_idx, JGA_date_idx, PFB_date_idx)
viable_Vinf_pluto(Launch_date_idx, JGA_date_idx, PFB_date_idx)
GA_vel_err_3d(Launch_date_idx, JGA_date_idx, PFB_date_idx)


out2.short_way_dv1_store(JGA_date_idx,PFB_date_idx)
out1.short_way_dv2_store(Launch_date_idx,JGA_date_idx)

%%
[r_earth ,v_earth] = MeeusEphemeris(Earth, Earth_dep_dates(Launch_date_idx), Sun);
[r_jup, v_jup] = MeeusEphemeris(Jupiter, Jupiter_arr_dates(JGA_date_idx), Sun);
[r_pluto, v_pluto] = MeeusEphemeris(Pluto, Pluto_arr_dates(PFB_date_idx), Sun); 

[~,v_jup_in] = lambert(r_earth, r_jup, ...
    (Jupiter_arr_dates(JGA_date_idx) ...
    - Earth_dep_dates(Launch_date_idx))*day2sec, 1, Sun);
[v_jup_out, ~] = lambert(r_jup, r_pluto, ...
    (Pluto_arr_dates(PFB_date_idx) ...
    - Jupiter_arr_dates(JGA_date_idx))*day2sec, 1, Sun);

[~, ~, ~, ~, rp] = BPlaneTarget(v_jup_in-v_jup, v_jup_out-v_jup, Jupiter.mu)
r_min = 2144760;
rp > r_min

%% 