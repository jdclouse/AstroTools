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
dur = final - begin

begin = zeros(1,32);
final = zeros(1,32);
ended = zeros(1,32);
dsite = nist;
for ii = 1:rows
    for jj = 1:cols
        if ~isnan(dsite.GPSdata_solar_max.topo_el(ii,jj)) && begin(jj) == 0 && final(jj) == 0
            begin(jj) = dsite.solar_max_hrofweek(ii);
        elseif ~isnan(dsite.GPSdata_solar_max.topo_el(ii,jj)) && begin(jj) ~= 0 && ~ended(jj)
            final(jj) = dsite.solar_max_hrofweek(ii);
        elseif isnan(dsite.GPSdata_solar_max.topo_el(ii,jj)) && begin(jj) ~= 0
            ended(jj) = 1;
        end
    end
end
dur = final - begin