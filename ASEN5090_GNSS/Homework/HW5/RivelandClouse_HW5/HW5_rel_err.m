%==========================================================================
%==========================================================================
% HW5_rel_err.m
%
% Author: Johnathan Clouse
%
% Relative error processing for HW5
%
%==========================================================================
%==========================================================================

%% Initialize
clear
close all

%% Read in GPS data in CSV format
% Ryan had the Delorme, with WAAS on
% John had the Garmin, WAAS off
[ryan_in, john_in] = import_gps_data('RyanInMat.csv', 'JohnInMat2.csv');
[ryan_out, john_out] = import_gps_data('RyanOut.csv', 'JohnOut.csv');

%% Convert errors...
% Compute arc lengths along the ellipsoid
% WGS84, http://earth-info.nga.mil/GandG/publications/tr8350.2/wgs84fin.pdf
a = 6378137.0; % m
b = 6356752.3142; % m
e = sqrt(a*a-b*b)/a;
lat = ryan_in.lat(1); % rad, just a representative latitude
lon2m = a*cos(lat)/sqrt(1-e*e*sin(lat)*sin(lat));
lat2m = a*(1-e*e)/(sqrt(1-e*e*sin(lat)*sin(lat))*(1-e*e*sin(lat)*sin(lat)));

err.lat = (john_in.lat-ryan_in.lat)*lat2m;
err.lon = (john_in.lon-ryan_in.lon)*lon2m;
err.alt = (john_in.alt-ryan_in.alt);
% err.rss = sqrt(err.lat.*err.lat+err.lon.*err.lon+err.alt.*err.alt);
err.rss = sqrt(err.lat.*err.lat+err.lon.*err.lon);
err.time = john_in.time;

%Waypoints
in_wps = ['04-OCT-14 20:28:09'
    '04-OCT-14 20:42:37'
    '04-OCT-14 20:51:42'
    '04-OCT-14 20:57:31'
    '04-OCT-14 21:05:10'
    '04-OCT-14 21:09:44'
    '04-OCT-14 21:17:12'
    '04-OCT-14 21:23:44'];

out_wps = ['04-OCT-14 21:23:44'
    '04-OCT-14 21:32:55'
    '04-OCT-14 21:37:04'
    '04-OCT-14 21:39:01'
    '04-OCT-14 21:43:37'
    '04-OCT-14 21:45:19'
    '04-OCT-14 21:50:58'
    '04-OCT-14 21:56:54'];

in_wps_times = zeros(7,1);
out_wps_times = zeros(7,1);
for ii = 1:8
    in_wps_times(ii) = datenum(in_wps(ii,:));
    out_wps_times(ii) = datenum(out_wps(ii,:));
end
in_wps_times(1)=err.time(1);

subplot(3,1,1)
lab_err_plots(err.time,err.lat,'UTC (hr)', 'Error (m)', 'Trip 1 Latitude Error', in_wps_times)

subplot(3,1,2)
lab_err_plots(err.time,err.lon,'UTC (hr)', 'Error (m)', 'Trip 1 Longitude Error', in_wps_times)

subplot(3,1,3)
lab_err_plots(err.time,err.rss,'UTC (hr)', 'Error (m)', 'Trip 1 RSS Error', in_wps_times)


errout.lat = (john_out.lat-ryan_out.lat)*lat2m;
errout.lon = (john_out.lon-ryan_out.lon)*lon2m;
errout.alt = (john_out.alt-ryan_out.alt);
% errout.rss = sqrt(errout.lat.*errout.lat+errout.lon.*errout.lon+errout.alt.*errout.alt);
errout.rss = sqrt(errout.lat.*errout.lat+errout.lon.*errout.lon);
errout.time = john_out.time;
out_wps_times(1) = errout.time(1);
out_wps_times = flip(out_wps_times); %so the waypoint plots show up correctly

figure
subplot(3,1,1)
lab_err_plots(errout.time,errout.lat,'UTC (hr)', 'Error (m)', 'Trip 2 Latitude Error', out_wps_times)

subplot(3,1,2)
lab_err_plots(errout.time,errout.lon,'UTC (hr)', 'Error (m)', 'Trip 2 Longitude Error', out_wps_times)

subplot(3,1,3)
lab_err_plots(errout.time,errout.rss,'UTC (hr)', 'Error (m)', 'Trip 2 RSS Error', out_wps_times)

figure
subplot(2,1,1)
lab_err_plots(err.time,err.alt,'UTC (hr)', 'Error (m)', 'Trip 1 Altitude Error', in_wps_times)
subplot(2,1,2)
lab_err_plots(errout.time,errout.alt,'UTC (hr)', 'Error (m)', 'Trip 2 Altitude Error', out_wps_times)