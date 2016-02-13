function [ fh ] = PorkchopPlot( DepartureDates, ArrivalDates, ...
    params)
%PorkchopPlot Create a porkchop plot
%   Detailed explanation goes here
fcnPrintQueue(mfilename('fullpath')) % Add this code to code app 

Sun = params.Sun;
planet1 = params.planet1;
planet2 = params.planet2;
c3_countours = params.c3_countours;
v_inf_countours = params.v_inf_countours;
TOF_countours = params.TOF_countours;
fig_dim = params.fig_dim;
day2sec = params.day2sec;

min_xfer_time = 100; %days
max_xfer_time = 500; %days

dep_elapsed_time = DepartureDates - DepartureDates(1);
arr_elapsed_time = ArrivalDates - ArrivalDates(1);

num_dep = length(DepartureDates);
num_arr = length(ArrivalDates);

short_way_dv1_store = NaN(num_dep, num_arr);
short_way_dv2_store = NaN(num_dep, num_arr);
long_way_dv1_store = NaN(num_dep, num_arr);
long_way_dv2_store = NaN(num_dep, num_arr);
TOF_store = zeros(num_dep, num_arr);

cx = 0;
% percent_complete = '';
for JD_depart = DepartureDates
    cx = cx + 1;
    TOF_store(cx,:) = ArrivalDates - DepartureDates(cx);
    cy = 0;    
%     fprintf(repmat('\b',size(percent_complete)))
%     percent_complete = [num2str(floor(cx/num_dep)) '%%'];
%     fprintf(percent_complete);
    for JD_arrive = ArrivalDates
        cy = cy + 1;
%         fprintf([num2str(cx) ', ' num2str(cy) '\n']);
        if JD_arrive-JD_depart < min_xfer_time
            continue % Don't bother with really short trajectories
        elseif JD_arrive-JD_depart > max_xfer_time
            continue % Don't bother with really long ones
        end
        [r_p1, v_p1] = MeeusEphemeris(planet1, JD_depart, Sun);
        [r_p2, v_p2] = MeeusEphemeris(planet2, JD_arrive, Sun);
        % Short way (Type I)
        [v1, v2] = lambert(r_p1, r_p2, (JD_arrive-JD_depart)*day2sec, ...
            1, Sun);
        short_way_dv1_store(cx, cy) = norm(v1 - v_p1);
        short_way_dv2_store(cx, cy) = norm(v2 - v_p2);
        % Long way (Type II)
        [v1, v2] = lambert(r_p1, r_p2, (JD_arrive-JD_depart)*day2sec, ...
            -1, Sun);
        long_way_dv1_store(cx, cy) = norm(v1 - v_p1);
        long_way_dv2_store(cx, cy) = norm(v2 - v_p2);
    end
end

sw_c3_store = short_way_dv1_store.*short_way_dv1_store;
lw_c3_store = long_way_dv1_store.*long_way_dv1_store;

fh = figure('Position', fig_dim);
hold on
[cs1, h1] = contour(dep_elapsed_time, arr_elapsed_time, sw_c3_store', ...
    c3_countours, 'r');
clabel(cs1, h1)
[cs2, h2] = contour(dep_elapsed_time, arr_elapsed_time, lw_c3_store', ...
    c3_countours, 'r');
clabel(cs2, h2)
[cs3, h3] = contour(dep_elapsed_time, arr_elapsed_time, short_way_dv2_store', ...
    v_inf_countours, 'b');
clabel(cs3, h3)
[cs4, h4] = contour(dep_elapsed_time, arr_elapsed_time, long_way_dv2_store', ...
    v_inf_countours, 'b');
clabel(cs4, h4)
[cs5, h5] = contour(dep_elapsed_time, arr_elapsed_time, TOF_store', ...
    TOF_countours, 'k');
clabel(cs5, h5)

legend([h1, h3, h5], {'C_3 (km^2/s^2)','Arrival V_{\infty} (km/s)',...
    'Time of Flight (days)'}, 'Location', 'NorthWest')
ax = gca;
% Warning, this may be weird for more precise applications. have to think
% about it.
jd_dn_delta = floor(juliandate(date)) - datenum(date);
dep_date_dn = DepartureDates(1) - jd_dn_delta;
arr_date_dn = ArrivalDates(1) - jd_dn_delta;

depDates = get(ax, 'XTick');
set(ax, 'XTickLabel', datestr(dep_date_dn+depDates,1));
% set(ax, 'XTickLabelRotation', 45);
arrDates = get(ax, 'YTick');
set(ax, 'YTickLabel', datestr(arr_date_dn+arrDates,1));
% set(ax, 'YTickLabelRotation', 45);

xlabel('Launch Date')
ylabel('Arrival Date')

fprintf('Plot generated.\n')
end
