%% Ionosphere Error Project Script
%% Initialize
%==========================================================================
% Simplified driver for GPS Visbility Codes
%
% Based on ASEN5090_GPSvisibility_script given to class by P. Axelrad 9/12
% Based on GPSVisibility_GUI by Ben K. Bradley and calling functions used
% in that GUI.
%
%==========================================================================

clearvars -except function_list pub_opt ntus nist solar_min solar_max
close all

if ~exist('ntus','var') || ~isfield(ntus,'obs_solar_min')
    ntus.obs_solar_min = read_rinex_obs4('ntus0790.09o');
end
if ~exist('ntus','var') || ~isfield(ntus,'obs_solar_max')
%     ntus.obs_solar_max = read_rinex_obs4('ntus0790.13o');
    ntus.obs_solar_max = read_rinex_obs4('ntus0790.14o');
end
if ~exist('nist','var') || ~isfield(nist,'obs_solar_min')
    nist.obs_solar_min = read_rinex_obs4('nist0790.09o');
end
if ~exist('nist','var') || ~isfield(nist,'obs_solar_max')
%     nist.obs_solar_max = read_rinex_obs4('nist0790.13o');
    nist.obs_solar_max = read_rinex_obs4('nist0790.14o');
end

%Site info
%SOURCE Singapore
ntus.latgd = 1.204488;  % latitude, deg
ntus.lon   = 103.404764;  % longitude, deg
ntus.alt   = 79.0;  % altitude, m
ntus.C1 = 6;
ntus.P2 = 7;
ntus.L1 = 4;
ntus.L2 = 5;
ntus.description = 'Singapore';

% NIST, Boulder, CO
% ftp://ftp.igs.org/pub/station/log/nist_20080725.log 
% Accessed 2013/10/26
nist.latgd = 39.594224;  % latitude, deg
nist.lon   = -105.154538;  % longitude, deg
nist.alt   = 1648.488;  % altitude, m
nist.C1 = 4;
nist.P2 = 8;
nist.L1 = 5;
nist.L2 = 9;
nist.description = 'NIST';

% Enter the time of interest
solar_min.UTC = [2009 3 20 0 0 0];
solar_min.idx = 1;
solar_min.description = 'Solar Min';
% solar_max.UTC = [2013 3 20 0 0 0];
solar_max.UTC = [2014 3 19 23 59 59]; % Leap-seconds = 16 now
solar_max.idx = 2;
solar_max.description = 'Solar Max';

% 2009
solar_min.a = [1.1694e-8 5.9461e-9 -3.6169e-07 4.1277e-07];
solar_min.b = [1.0511e05 -4.1263e04 -1.3974e06  5.4206e06]';
% 2013
solar_max.a = [2.8197e-08  2.4052e-08 -2.7484e-07 -4.7981e-07 ];
solar_max.b = [1.6588e05 -1.8886e05 -3.3702e06  1.3382e07]';
% 2014
solar_max.a = [3.5404e-08  1.1407e-08 -1.3468e-07 -4.3045e-08 ];
solar_max.b = [1.7226e05  1.3073e05 -8.6663e05 -2.1569e06 ]';

solar_max.ionex = read_ionex_obs('jplg0790.14i');
solar_max.ionex.epochs(end) = solar_max.ionex.epochs(end-1) + 7200;
solar_min.ionex = read_ionex_obs('jplg0790.09i');
solar_min.ionex.epochs(end) = solar_min.ionex.epochs(end-1) + 7200;

% Grab duration and time step =============================================
durationhrs = 24;

dt_sec = 15;

% Input Antenna Location ===================================================
% Make antenna pointed straight up
ant_enu = [0 0 1];

% Set minimum and maximum mask angles
mask_min = 0;  % deg
mask_max = 90; % deg

if ~isfield(ntus, 'GPSdata_solar_min')
    solar_min.GPSvec = utc2gpsvec(solar_min.UTC);  % This will adjust for leap second offset
    [WN2, TOW, WN1] = gpsvec2gpstow(solar_min.GPSvec) ; %WN1 is full week, WN2 is mod1024
    check = gpstow2gpsvec(WN2,TOW,2); % See if it converts back correctly, 2 for WN2
    % Construct YUMA almanac file name since this is default setting
    solar_min.navfilename = generate_GPSyuma_name(solar_min.GPSvec);
    [solar_min.navfilename,statusflag] = download_GPSyuma(solar_min.GPSvec);
    
    [ntus.solar_min_time_wntow,ntus.GPSdata_solar_min] = ...
        ASEN5090_GPSvis(solar_min.navfilename, 1, solar_min.GPSvec,...
        durationhrs, dt_sec, ntus.latgd, ntus.lon, ntus.alt,...
        mask_min, mask_max, mask_min, ant_enu, 0, []);
    ntus.solar_min_hrofweek = ntus.solar_min_time_wntow(:,2)/3600;
end

if ~isfield(ntus, 'GPSdata_solar_max')
    solar_max.GPSvec = utc2gpsvec(solar_max.UTC);  % This will adjust for leap second offset
    [WN2, TOW, WN1] = gpsvec2gpstow(solar_max.GPSvec) ; %WN1 is full week, WN2 is mod1024
    check = gpstow2gpsvec(WN2,TOW,2); % See if it converts back correctly, 2 for WN2
    % Construct YUMA almanac file name since this is default setting
    solar_max.navfilename = generate_GPSyuma_name(solar_max.GPSvec);
    [solar_max.navfilename,statusflag] = download_GPSyuma(solar_max.GPSvec);
    
    [ntus.solar_max_time_wntow,ntus.GPSdata_solar_max] = ...
        ASEN5090_GPSvis(solar_max.navfilename, 1, solar_max.GPSvec,...
        durationhrs, dt_sec, ntus.latgd, ntus.lon, ntus.alt,...
        mask_min, mask_max, mask_min, ant_enu, 0, []);
    ntus.solar_max_hrofweek = ntus.solar_max_time_wntow(:,2)/3600;
end

if ~isfield(nist, 'GPSdata_solar_min')
    solar_min.GPSvec = utc2gpsvec(solar_min.UTC);  % This will adjust for leap second offset
    [WN2, TOW, WN1] = gpsvec2gpstow(solar_min.GPSvec) ; %WN1 is full week, WN2 is mod1024
    check = gpstow2gpsvec(WN2,TOW,2); % See if it converts back correctly, 2 for WN2
    % Construct YUMA almanac file name since this is default setting
    solar_min.navfilename = generate_GPSyuma_name(solar_min.GPSvec);
    [solar_min.navfilename,statusflag] = download_GPSyuma(solar_min.GPSvec);
    
    [nist.solar_min_time_wntow,nist.GPSdata_solar_min] = ...
        ASEN5090_GPSvis(solar_min.navfilename, 1, solar_min.GPSvec,...
        durationhrs, dt_sec, nist.latgd, nist.lon, nist.alt,...
        mask_min, mask_max, mask_min, ant_enu, 0, []);
    nist.solar_min_hrofweek = nist.solar_min_time_wntow(:,2)/3600;
end

if ~isfield(nist, 'GPSdata_solar_max')
    solar_max.GPSvec = utc2gpsvec(solar_max.UTC);  % This will adjust for leap second offset
    [WN2, TOW, WN1] = gpsvec2gpstow(solar_max.GPSvec) ; %WN1 is full week, WN2 is mod1024
    check = gpstow2gpsvec(WN2,TOW,2); % See if it converts back correctly, 2 for WN2
    % Construct YUMA almanac file name since this is default setting
    solar_max.navfilename = generate_GPSyuma_name(solar_max.GPSvec);
    [solar_max.navfilename,statusflag] = download_GPSyuma(solar_max.GPSvec);
    
    [nist.solar_max_time_wntow,nist.GPSdata_solar_max] = ...
        ASEN5090_GPSvis(solar_max.navfilename, 1, solar_max.GPSvec,...
        durationhrs, dt_sec, nist.latgd, nist.lon, nist.alt,...
        mask_min, mask_max, mask_min, ant_enu, 0, []);
    nist.solar_max_hrofweek = nist.solar_max_time_wntow(:,2)/3600;
end

% [time_wntow,GPSdata] = ASEN5090_GPSvis(navfilename, 1, GPSvec,...
%     durationhrs, dt_sec, latgd, lon, alt,...
%     mask_min, mask_max, mask_min, ant_enu, 0, []);
% hrofweek = time_wntow(:,2)/3600;

% Plot results ============================================================
% =========================================================================


% Number of Satellites Visible ---------------------------
[ax2] = plot_GPSnumsats(ntus.solar_min_hrofweek,ntus.GPSdata_solar_min.ant_numsats);

title(ax2, 'Number of visible satellites')

% Topocentric: AzEl Plot --------------------------------------------------
[rows,cols] = size(ntus.GPSdata_solar_min.topo_el);
az_vec     = reshape(ntus.GPSdata_solar_min.topo_az,rows*cols,1);
el_vec     = reshape(ntus.GPSdata_solar_min.topo_el,rows*cols,1);

% ntus.GPSdata_solar_min.prn = repmat([1:32],rows,1);

% Find the longest times
begin = zeros(1,32);
final = zeros(1,32);
ended = zeros(1,32);
dsite = nist;
for ii = 1:rows
    for jj = 1:cols
        if ~isnan(dsite.GPSdata_solar_min.topo_el(ii,jj)) && begin(jj) == 0 && final(jj) == 0
            begin(jj) = dsite.solar_min_hrofweek(ii);
        elseif ~isnan(dsite.GPSdata_solar_min.topo_el(ii,jj)) && begin(jj) ~= 0 && ~ended(jj)
            final(jj) = dsite.solar_min_hrofweek(ii);
        elseif isnan(dsite.GPSdata_solar_min.topo_el(ii,jj)) && begin(jj) ~= 0
            ended(jj) = 1;
        end
    end
end
dur = final - begin;

prns_nist = [14 18 29];
prns_ntus = [12 14 22];    

% plot the sats visible to antenna only%fig3 = figure; ax3 = axes;
figure
% plotAzEl(GPSdata.topo_az',GPSdata.topo_el',GPSdata.prn')
% plotAzEl(ntus.GPSdata_solar_min.topo_az(:,12)',ntus.GPSdata_solar_min.topo_el(:,12)',ntus.GPSdata_solar_min.prn(:,12)')
% hold on
% plotAzEl(ntus.GPSdata_solar_min.topo_az(:,14)',ntus.GPSdata_solar_min.topo_el(:,14)',ntus.GPSdata_solar_min.prn(:,14)')
% plotAzEl(ntus.GPSdata_solar_min.topo_az(:,32)',ntus.GPSdata_solar_min.topo_el(:,32)',ntus.GPSdata_solar_min.prn(:,32)')
plotAzEl(ntus.GPSdata_solar_min.topo_az(:,12)',ntus.GPSdata_solar_min.topo_el(:,12)')
hold on
plotAzEl(ntus.GPSdata_solar_min.topo_az(:,14)',ntus.GPSdata_solar_min.topo_el(:,14)')
plotAzEl(ntus.GPSdata_solar_min.topo_az(:,32)',ntus.GPSdata_solar_min.topo_el(:,22)')

% Antenna-centric: Elevation Plot -----------------------------------------
msize = [rows cols];

el_plots(prns_ntus, ntus.GPSdata_solar_min.topo_el, ntus.solar_min_hrofweek, ...
    'Solar Min', 'Singapore', msize)
el_plots(prns_ntus, ntus.GPSdata_solar_max.topo_el, ntus.solar_max_hrofweek, ...
    'Solar Max', 'Singapore', msize)
el_plots(prns_nist, nist.GPSdata_solar_min.topo_el, nist.solar_min_hrofweek, ...
    'Solar Min', 'NIST', msize)
el_plots(prns_nist, nist.GPSdata_solar_max.topo_el, nist.solar_max_hrofweek, ...
    'Solar Max', 'NIST', msize)

%%
figure
hold on
for ii = 1:32
    % For a given PRN...
    prn_indices_el = ~isnan(ntus.GPSdata_solar_min.topo_el(:,ii));
    temp_el_times = ntus.solar_min_time_wntow(prn_indices_el,2);
    temp_el = ntus.GPSdata_solar_min.topo_el(prn_indices_el,ii);
    prn_indices_data = ntus.obs_solar_min.data(:,3) == ii;
    temp_data = ntus.obs_solar_min.data(prn_indices_data,:);
    % For all common measurement/elevation times...
    [C, ai, bi]=intersect(temp_data(:,2),temp_el_times);
    % Compute Iono delay...
    if ~isempty(temp_el)
    [ Idelay, Iz ] = df_iono_delay( temp_data(ai,6), ...
        temp_data(ai,7), temp_el(bi)*(pi/180));
    subplot(2,1,1)
    hold on
    plot(temp_el_times(bi)/3600, Idelay)
    subplot(2,1,2)
    hold on
    plot(temp_el(bi), Idelay)
    end        
end

figure
hold on
for ii = 1:32
    % For a given PRN...
    prn_indices_el = ~isnan(ntus.GPSdata_solar_max.topo_el(:,ii));
    temp_el_times = ntus.solar_max_time_wntow(prn_indices_el,2);
    temp_el = ntus.GPSdata_solar_max.topo_el(prn_indices_el,ii);
    prn_indices_data = ntus.obs_solar_max.data(:,3) == ii;
    temp_data = ntus.obs_solar_max.data(prn_indices_data,:);
    % For all common measurement/elevation times...
    [C, ai, bi]=intersect(temp_data(:,2),temp_el_times);
    % Compute Iono delay...
    if ~isempty(temp_el)
    [ Idelay, Iz ] = df_iono_delay( temp_data(ai,6), ...
        temp_data(ai,7), temp_el(bi)*(pi/180));
    subplot(2,1,1)
    hold on
    plot(temp_el_times(bi)/3600, Idelay)
    subplot(2,1,2)
    hold on
    plot(temp_el(bi), Idelay)
    end        
end

%%
% trim some times from the NIST @ solar min to get better data viz
nist.GPSdata_solar_min.topo_az((nist.solar_min_hrofweek - nist.solar_min_hrofweek(1))<6,14) = NaN;
nist.GPSdata_solar_min.topo_el((nist.solar_min_hrofweek - nist.solar_min_hrofweek(1))<6,14) = NaN;
nist.GPSdata_solar_min.topo_az((nist.solar_min_hrofweek - nist.solar_min_hrofweek(1))<6,29) = NaN;
nist.GPSdata_solar_min.topo_el((nist.solar_min_hrofweek - nist.solar_min_hrofweek(1))<6,29) = NaN;
plot_sat_data(prns_ntus, ntus, ntus.GPSdata_solar_min, ntus.solar_min_time_wntow, ntus.obs_solar_min, solar_min)
plot_sat_data(prns_ntus, ntus, ntus.GPSdata_solar_max, ntus.solar_max_time_wntow, ntus.obs_solar_max, solar_max)
plot_sat_data(prns_nist, nist, nist.GPSdata_solar_min, nist.solar_min_time_wntow, nist.obs_solar_min, solar_min)
plot_sat_data(prns_nist, nist, nist.GPSdata_solar_max, nist.solar_max_time_wntow, nist.obs_solar_max, solar_max)
