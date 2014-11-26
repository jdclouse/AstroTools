function TEC = iono_interp(lat, lon, t, ionex)
%TODO: error checking for polar conditions.
% Uses format definition from 
%   IONEX: The IONosphere Map EXchange Format Version 1 (1998)
%   Schaer, Gurtner, Feltens

% find the indices of lat/lon in the ionex info between which is your
% location
% Assume maps are given North to South, West to East
lat0 = 0;
lat1 = 0;
dlat = 0;
lat0_idx = -1;
lat1_idx = -1;
for ii = 1:length(ionex.lats)
    if ionex.lats(ii) <= lat
        lat0_idx = ii;
        lat1_idx = ii-1;
        lat0 = ionex.lats(ii); % less than input
        lat1 = ionex.lats(ii-1); % greater than input
        dlat = lat1-lat0;
        break
    end
end
    
lon0 = 0;
lon1 = 0;
dlon = 0;
lon0_idx = -1;
lon1_idx = -1;
for ii = 1:length(ionex.lons)
    if ionex.lons(ii) >= lon
        lon1_idx = ii;
        lon0_idx = ii-1;
        lon1 = ionex.lons(ii); % less than input
        lon0 = ionex.lons(ii-1); % greater than input
        dlon = lon1-lon0;
        break
    end
end    

t0 = 0;
t1 = 0;
dt = 0;
t0_idx = -1;
t1_idx = -1;
for ii = 1:length(ionex.epochs)
    if t <= ionex.epochs(ii)
        if ii == 1
            t0_idx = ii;
            t1_idx = ii;
            t0 = ionex.epochs(ii); % less than input
            t1 = ionex.epochs(ii); % less than input
        else
            t0_idx = ii-1;
            t1_idx = ii;
            t0 = ionex.epochs(ii-1); % less than input
            t1 = ionex.epochs(ii); % less than input
        end
        dt = t1-t0;
        break
    end
end

% From Schaer et al.
% Bivariate interpolation for TEC on the map at epochs.
% Linear interpolation between the epochs.
q = (lat-lat0)/dlat;
p = (lon-lon0)/dlon;
tdiff = t-t0;

E00 = ionex.maps(t0_idx,lat0_idx,lon0_idx);
E10 = ionex.maps(t0_idx,lat0_idx,lon1_idx);
E01 = ionex.maps(t0_idx,lat1_idx,lon0_idx);
E11 = ionex.maps(t0_idx,lat1_idx,lon1_idx);
TEC0 = (1-p)*(1-q)*E00 + p*(1-q)*E10 + q*(1-p)*E01 + p*q*E11;

E00 = ionex.maps(t1_idx,lat0_idx,lon0_idx);
E10 = ionex.maps(t1_idx,lat0_idx,lon1_idx);
E01 = ionex.maps(t1_idx,lat1_idx,lon0_idx);
E11 = ionex.maps(t1_idx,lat1_idx,lon1_idx);
TEC1 = (1-p)*(1-q)*E00 + p*(1-q)*E10 + q*(1-p)*E01 + p*q*E11;

m = (TEC1-TEC0)/dt;
if isnan(m)
    m = 0;
end
TEC = TEC0 + tdiff*m;