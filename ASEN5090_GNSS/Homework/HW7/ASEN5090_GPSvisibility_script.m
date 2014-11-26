%==========================================================================
% Simplified driver for GPS Visbility Codes
%
% Based on GPSVisibility_GUI by Ben K. Bradley and calling functions used
% in that GUI.
% P. Axelrad 9/12
%

%==========================================================================
%==========================================================================

% clear, close all

% Enter the time of interest
% UTC = [2014 9 27 19 0 0];
UTC = [2009 3 20 4 0 0];
% UTC = [2014 9 20 19 0 0];
GPSvec = utc2gpsvec(UTC);  % This will adjust for leap second offset
[WN2, TOW, WN1] = gpsvec2gpstow(GPSvec) ; %WN1 is full week, WN2 is mod1024
check = gpstow2gpsvec(WN2,TOW,2); % See if it converts back correctly, 2 for WN2


% Construct YUMA almanac file name since this is default setting
navfilename = generate_GPSyuma_name(GPSvec);
[navfilename,statusflag] = download_GPSyuma(GPSvec);


% Grab duration and time step =============================================
durationhrs = 5;

dt_sec = 600;


% Input Antenna Location ===================================================

%SOURCE
latgd = 1.204488;  % latitude, deg
lon   = 103.404764;  % longitude, deg
alt   = 79.0;  % altitude, m

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
[ax2] = plot_GPSnumsats(hrofweek,GPSdata.ant_numsats);

title(ax2, 'Number of visible satellites')

% Topocentric: AzEl Plot --------------------------------------------------
[rows,cols] = size(GPSdata.topo_el);

first = 1;
last = 3;
az_vec     = reshape(GPSdata.topo_az,rows*cols,1);
el_vec     = reshape(GPSdata.topo_el,rows*cols,1);

GPSdata.prn = repmat([1:32],rows,1);

% Find the longest times
begin = zeros(1,32);
final = zeros(1,32);
ended = zeros(1,32);
dur = zeros(1,32);
for ii = 1:rows
    for jj = 1:cols
        if ~isnan(GPSdata.topo_el(ii,jj)) && begin(jj) == 0 && final(jj) == 0
            begin(jj) = hrofweek(ii);
        elseif ~isnan(GPSdata.topo_el(ii,jj)) && begin(jj) ~= 0 && ~ended(jj)
            final(jj) = hrofweek(ii);
        elseif isnan(GPSdata.topo_el(ii,jj)) && begin(jj) ~= 0
            ended(jj) = 1;
        end
    end
end
dur = final - begin;
    

% plot the sats visible to antenna only%fig3 = figure; ax3 = axes;
figure
plotAzEl(GPSdata.topo_az',GPSdata.topo_el',GPSdata.prn')
% plotAzEl(GPSdata.topo_az(:,12)',GPSdata.topo_el(:,12)',GPSdata.prn(:,12)')
% hold on
% plotAzEl(GPSdata.topo_az(:,14)',GPSdata.topo_el(:,14)',GPSdata.prn(:,14)')
% plotAzEl(GPSdata.topo_az(:,32)',GPSdata.topo_el(:,32)',GPSdata.prn(:,32)')

prn_i_want = [12 14 32];
masked_el = nan(rows,cols);
for ii = 1:length(prn_i_want)
    masked_el(:,prn_i_want(ii)) = GPSdata.topo_el(:,prn_i_want(ii));
end
masked_el_vec     = reshape(masked_el,rows*cols,1);

% Antenna-centric: Elevation Plot -----------------------------------------

time_mat = repmat(hrofweek,1,cols);
time_vec = reshape(time_mat,rows*cols,1);

fig4 = figure; ax4 = axes;
% plot(ax4,time_vec,el_vec,'ob','markerfacecolor','b','markersize',4);
plot(ax4,time_vec,masked_el_vec,'ob','markerfacecolor','b','markersize',4);

ylabel('Elevation (deg)');
xlabel('Time (hr)');

grid(ax4,'on');


title(ax4,{'\fontsize{11}\bfElevation Angle';'\bfof Satellites Seen by Antenna'});

%% CLOUSE: Load 9/27 data and subtract to find the satellite az/el diff

% GPSdata_new=GPSdata;
% load('GPS927.mat')
% %Display the az/el diffs
% GPSdata_new.prn(1,GPSdata_new.topo_el(1,:)>0)
% GPSdata_new.topo_az(1,GPSdata_new.topo_el(1,:)>0)-GPSdata.topo_az(1,GPSdata_new.topo_el(1,:)>0)
% GPSdata_new.topo_el(1,GPSdata_new.topo_el(1,:)>0)-GPSdata.topo_el(1,GPSdata_new.topo_el(1,:)>0)
