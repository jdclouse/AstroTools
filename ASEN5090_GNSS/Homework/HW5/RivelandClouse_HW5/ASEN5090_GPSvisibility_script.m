%==========================================================================
% Simplified driver for GPS Visbility Codes
%
% Based on GPSVisibility_GUI by Ben K. Bradley and calling functions used
% in that GUI.
% P. Axelrad 9/12
%

%==========================================================================
%==========================================================================

clear, close all

% Enter the time of interest

%Local computation

%Wpt1
%UTC = [2014 10 4 20 30 26]
% %Wpt2
%UTC = [2014 10 4 20 42 35]
% %Wpt3
 %UTC = [2014 10 4 20 52 34]
% %Wpt4
% UTC = [2014 10 4 20 57 07]
%Wpt5
 %UTC = [2014 10 4 21 04 55]
%Wpt6
 %UTC = [2014 10 4 21 09 40]
%Wpt7
%UTC = [2014 10 4 21 19 07]
%Wpt8
UTC = [2014 10 4 21 23 29]

GPSvec = utc2gpsvec(UTC);  % This will adjust for leap second offset
[WN2, TOW, WN1] = gpsvec2gpstow(GPSvec) ; %WN1 is full week, WN2 is mod1024
check = gpstow2gpsvec(WN2,TOW,2); % See if it converts back correctly, 2 for WN2

% Construct YUMA almanac file name since this is default setting
navfilename = generate_GPSyuma_name(GPSvec);
[navfilename,statusflag] = download_GPSyuma(GPSvec);

% Grab duration and time step =============================================
durationhrs = 24;
dt_sec = 600;

% Input Antenna Location ===================================================
% My house location
% latgd = -40.966066;  % latitude, deg
% lon   = -104.720037;  % longitude, deg
% alt   = 2142.744;  % altitude, m

% %Wpt1
% latgd = 39.427197695;  % latitude, deg
% lon   = -104.902699351;  % longitude, deg
% alt   = 1836.196167;  % altitude, m

% %Wpt2
% latgd = 39.428788066;  % latitude, deg
% lon   = -104.969208837;  % longitude, deg
% alt   = 1766.585815;  % altitude, m
% 
% %Wpt3
% latgd = 39.386439085;  % latitude, deg
% lon   = -105.022118211;  % longitude, deg
% alt   = 1951.067993;  % altitude, m
% % 
%Wpt4
% latgd = 39.383223891;  % latitude, deg
% lon   = -105.031314731;  % longitude, deg
% alt   = 2039.067139;  % altitude, m
% 
% %Wpt5
% latgd = 39.385950208;  % latitude, deg
% lon   = -105.075039029;  % longitude, deg
% alt   = 2217.168213;  % altitude, m
% 
%Wpt6
% latgd = 39.383643746;  % latitude, deg
% lon   = -105.075039029;  % longitude, deg
% alt   = 2217.168213;  % altitude, m
% 
% %Wpt7
% latgd = 39.372378588;  % latitude, deg
% lon   = -105.106248617;  % longitude, deg
% alt   = 2225.073975;  % altitude, m
% 
% %Wpt8
latgd = 39.349440455;  % latitude, deg
lon   = -105.120900393;  % longitude, deg
alt   = 2163.671875;  % altitude, m

% Make antenna pointed straight up
ant_enu = [0 0 1];

% Set minimum and maximum mask angles
mask_min = 0;  % deg
mask_max = 90; % deg

[time_wntow,GPSdata] = ASEN5090_GPSvis(navfilename, 1, GPSvec,...
    durationhrs, dt_sec, latgd, lon, alt,...
    mask_min, mask_max, mask_min, ant_enu, 0, []);
hrofweek = time_wntow(:,2)/3600;

% Plot results ============================================================
% =========================================================================


% Number of Satellites Visible ---------------------------
% [ax2] = plot_GPSnumsats(hrofweek,GPSdata.ant_numsats);
% 
% title(ax2, 'Number of visible satellites')

% Topocentric: AzEl Plot --------------------------------------------------
[rows,cols] = size(GPSdata.topo_el);

az_vec     = reshape(GPSdata.topo_az,rows*cols,1);
el_vec     = reshape(GPSdata.topo_el,rows*cols,1);
GPSdata.prn = repmat([1:32],rows,1);

%plot the prns at beginning of time period, this polar plot shows az and el
% with labeled prns

%Compute WAAS vehicle az and el
Radius = 42164000;
AzOrbit = 252.7 * pi/180;
ElOrbit = 0.0;

Ecef = zeros(1,3);
Ecef(1) = Radius * cos(AzOrbit);
Ecef(2) = Radius * sin(AzOrbit);
Ecef(3) = 0.0;
% Compute ECEF Position of Antenna Location ===============================
r_site = lla2ecef(latgd, lon, alt*0.001); 

r_site = r_site' * 1000;  % [x y z] meters
[azWAAS, elWAAS, RangeWAAS] = ASEN5090_ecef2azelrange(Ecef,r_site,latgd,lon);

figure 
az = [GPSdata.topo_az(1,:),azWAAS];
el = [GPSdata.topo_el(1,:),elWAAS];
prn = [GPSdata.prn(1,:),138];
plotAzEl(az,el,prn,'*k')

[az,el] = ComputeElMask(8);
prns = zeros(1,360);
plotAzEl(az,el,prns,'*b');
% 
% title('Azimuth vs. Elevation, Colorado Springs, CO \newline                Sep 26, 2004 0300 UTC', 'fontsize',16);

% plot the sats visible to antenna only%fig3 = figure; ax3 = axes;
% figure
% plotAzEl(GPSdata.topo_az',GPSdata.topo_el',zeros(rows,cols))

%title('Azimuth vs. Elevation, Waypoint 1', 'fontsize',16);
 %title('Azimuth vs. Elevation, Waypoint 2', 'fontsize',16);
% title('Azimuth vs. Elevation, Waypoint 3', 'fontsize',16);
%title('Azimuth vs. Elevation, Waypoint 4', 'fontsize',16);
%title('Azimuth vs. Elevation, Waypoint 5', 'fontsize',16);
 %title('Azimuth vs. Elevation, Waypoint 6', 'fontsize',16);
  %title('Azimuth vs. Elevation, Waypoint 7', 'fontsize',16);
 title('Azimuth vs. Elevation, Waypoint 8', 'fontsize',16);

% figure
% plot(hrofweek,GPSdata.topo_el(:,5),hrofweek,GPSdata.topo_el(:,29),hrofweek,...
%     GPSdata.topo_el(:,17),hrofweek,GPSdata.topo_el(:,27));
% 
% ylabel('Elevation (deg)','fontsize',16);
% xlabel('Time of Week (hr)','fontsize',16);
% 
% title('Elevation vs. Time of Week', 'fontsize',16);
% legend('PRN 5','PRN 29','PRN 17','PRN 27');




%Antenna-centric: Elevation Plot -----------------------------------------

% time_mat = repmat(hrofweek,1,cols);
% time_vec = reshape(time_mat,rows*cols,1);
% 
% fig4 = figure; ax4 = axes;
% plot(ax4,time_vec,el_vec,'ob','markerfacecolor','b','markersize',4);
% 
% ylabel('Elevation (deg)');
% xlabel('Time (hr)');
% 
% grid(ax4,'on');
% 
% 
% title(ax4,{'\fontsize{11}\bfElevation Angle';'\bfof Satellites Seen by Antenna'});



