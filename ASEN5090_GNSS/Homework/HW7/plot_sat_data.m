function plot_sat_data(prn_i_want, site, sat_data, data_times, measurements,...
    solar_data)

TEC_IPP_plots = 0;
msz = 6;
sec_in_day = 86400;
counter = 1;
vs_el = figure;
subplot(1,3,counter);
vs_lt = figure;
subplot(1,3,counter);
ylabel('Ionospheric Delay (m)')
evt = figure;
subplot(3,1,counter);
azel = figure;
subplot(3,1,counter);
if TEC_IPP_plots
    llt = figure;
    subplot(3,1,counter);
    tec = figure;
    subplot(3,1,counter);
end
C1 = site.C1;
P2 = site.P2;
L1 = site.L1;
L2 = site.L2;
f_L1 = 1575.42e6; % MHz
f_L2 = 1227.60e6; % MHz
llh = [
    site.latgd*pi/180
    site.lon*pi/180
    site.alt];
for ii = prn_i_want
    % For a given PRN...
    prn_indices_el = ~isnan(sat_data.topo_el(:,ii));
    temp_el_times = data_times(prn_indices_el,2);
    temp_el = sat_data.topo_el(prn_indices_el,ii);
    temp_az = sat_data.topo_az(prn_indices_el,ii);
    prn_indices_data = measurements.data(:,3) == ii;
    temp_data = measurements.data(prn_indices_data,:);
    % For all common measurement/elevation times...
    [~, ai, bi]=intersect(temp_data(:,2),temp_el_times);
    
    % Local time
    cmn_az = temp_az(bi)*(pi/180);
    cmn_el = temp_el(bi)*(pi/180);
    cmn_utc_time = temp_el_times(bi); % TODO is this right?
    % Bound to a UTC day
%     while (max(cmn_time) > sec_in_day || min(cmn_time) < 0)
    while (max(cmn_utc_time) > sec_in_day) % assumes full-day's data, no more
%     ind = cmn_time < 0;
%     cmn_time = [cmn_time(~ind); cmn_time(ind)+sec_in_day];
        ind = cmn_utc_time >sec_in_day;
        cmn_utc_time = [cmn_utc_time(ind)-sec_in_day; cmn_utc_time(~ind)];
    end
%     lcl_time = (cmn_utc_time);
    lcl_time = zeros(length(cmn_utc_time),1);
    ipp_lat = zeros(length(cmn_utc_time),1);
    ipp_lon = zeros(length(cmn_utc_time),1);
    I_delay_map = zeros(length(cmn_utc_time),1);
    TEC = zeros(length(cmn_utc_time),1);
    for jj = 1:length(temp_el_times(bi))
        [lcl_time(jj), ipp_lat(jj), ipp_lon(jj)] = ...
            ipp_local_time( llh, cmn_az(jj), ...
            cmn_el(jj), cmn_utc_time(jj));
        
        TEC(jj) = iono_interp(ipp_lat(jj), ipp_lon(jj), cmn_utc_time(jj), ...
            solar_data.ionex);
        I_delay_map(jj) = ...
            40.3*...
            iono_interp(ipp_lat(jj), ipp_lon(jj), cmn_utc_time(jj), ...
            solar_data.ionex)/(f_L1*f_L1)*1e16;
    end
    
    % TEC and IPP Plots
    if TEC_IPP_plots
        figure(llt)
        subplot(3,1,counter)
        hold on
        plot(lcl_time/3600, ipp_lat*180/pi, 'r.')
        plot(lcl_time/3600, ipp_lon*180/pi, 'b.')
        figure(tec)
        subplot(3,1,counter)
        hold on
        plot(lcl_time/3600, TEC, 'b.')
    end
    
    % Compute Iono delay...
    [ Idelay, Iz ] = df_iono_delay( temp_data(ai,C1), ...
        temp_data(ai,P2), temp_el(bi)*(pi/180));
    figure(vs_el)
    subplot(1,3,counter)
    hold on
    plot(temp_el(bi), Idelay, '.', 'MarkerSize', msz)
    
    figure(vs_lt)
    subplot(1,3,counter)
    hold on
    plot(lcl_time/3600, Iz, '.', 'MarkerSize', msz)
    
    % Carrier divergence
    I_code_carrier_div = codecarrier_iono_divergence(...
        temp_data(ai,C1),...
        temp_data(ai,L1));
    % Adjust the c-c div to account for different passes
    for pp = 2:length(I_code_carrier_div)
        if lcl_time(pp)-lcl_time(pp-1) > 3600
            I_code_carrier_div(pp:end) = I_code_carrier_div(pp:end)...
                -I_code_carrier_div(pp);
        end
    end
    figure(vs_el)
    plot(temp_el(bi), I_code_carrier_div, 'r.', 'MarkerSize', msz);
    figure(vs_lt)
    plot(lcl_time/3600, I_code_carrier_div./iono_obliq_factor(cmn_el), 'r.'...
        , 'MarkerSize', msz)
    
    % Klobuchar
    [ Idelay, Iz, t_local ] = ...
        klobuchar( temp_az(bi)*(pi/180), temp_el(bi)*(pi/180), ...
        site.latgd*(pi/180), site.lon*(pi/180), ...
        cmn_utc_time, solar_data.a, solar_data.b );
    figure(vs_el)
    plot(temp_el(bi), Idelay, 'g.', 'MarkerSize', msz)
    figure(vs_lt)
    plot(t_local/3600, Iz, 'g.', 'MarkerSize', msz)

    % TEC map
    figure(vs_el)
    plot(temp_el(bi), I_delay_map.*iono_obliq_factor(cmn_el), 'c.', 'MarkerSize', msz)
    figure(vs_lt)
    plot(t_local/3600, I_delay_map, 'c.', 'MarkerSize', msz)
    
    % Plot elevations vs time where there was RINEX data
    figure(evt);
    subplot(3,1,counter);
    plot(lcl_time/3600, cmn_el*180/pi, '.')
    
    figure(azel);
    subplot(3,1,counter);
    plotAzEl(temp_az(bi),temp_el(bi))
    
    counter = counter + 1;
end

when = sprintf('\\bfDuring %s', solar_data.description);
where = sprintf('\\bfat %s', site.description);
prn_title = @(x) sprintf('PRN %d',prn_i_want(x));
sp_title = @(x) title(prn_title(x));

figure(vs_el);
subplot(1,3,2);
title({'\fontsize{11}\bfIono Delay vs. Elevation';when;where});
figure(vs_lt);
subplot(1,3,2);
title({'\fontsize{11}\bfIono Zenith Delay vs. IPP Local Time';when;where});
figure(evt);
subplot(3,1,1);
title({'\fontsize{11}\bfSat Elevation vs. IPP Local Time';when;where});

for ii = 1:3
    figure(vs_el);
    subplot(1,3,ii);
    xlabel(sprintf('%s Elevation (deg)', prn_title(ii)))
    figure(vs_lt);
    subplot(1,3,ii);
%     subplot(3,1,ii);
    xlabel(sprintf('%s IPP Local Time (hr)', prn_title(ii)))
    figure(evt);
    subplot(3,1,ii);
    xlabel(sprintf('%s IPP Local Time (hr)', prn_title(ii)))
    figure(azel);
    subplot(3,1,ii);
    title(sprintf('%s Az-El, %s %s', prn_title(ii), site.description, solar_data.description))
end
figure(vs_el);
legend('DF Code Delay','Code-Carrier Divergence','Klobuchar','TEC Map',...
    'Location', 'Best')
figure(vs_lt);
legend('DF Code Delay','Code-Carrier Divergence','Klobuchar','TEC Map',...
    'Location', 'Best')