function topo = ecef2topo( r_ecef, lat, lon, alt)
%ecef2topo return topocentric parameters from ECEF position, and the 
% origin's geocentric latitude, longitude, and altitude.
% Topo params are azimuth, elevation, and range (km)
% angle units in radians
fcnPrintQueue(mfilename('fullpath')) % Add this code to code app 

% Get some useful constants
CelestialConstants;

% ECEF of tracking station/topo origin
O_ECEF = latlonalt2ECEF(lat, lon, alt);

r_topo = norm(r_ecef - O_ECEF);

% convert both r_ecef and O_ECEF to East-North-Up coords (Misra/Enge)
% Ideally this should be done with geodetic latitude... FIXME later
ECEF2ENU_DCM = Euler2DCM('31',[lon+pi/2, pi/2-lat]);
r_enu = ECEF2ENU_DCM*(r_ecef - O_ECEF);
el = asin(r_enu(3)/norm(r_enu));
az = atan2(r_enu(1),r_enu(2));
if az < 0
    az = az + 2*pi;
end

topo = [az; el; r_topo];