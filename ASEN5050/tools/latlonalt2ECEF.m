function r_ECEF = latlonalt2ECEF( lat, lon, alt)
%latlonalt2ECEF return ECEF coords from geocentric latitude, longitude, 
% and altitude.
% units in radians
fcnPrintQueue(mfilename('fullpath')) % Add this code to code app 

% Get some useful constants
CelestialConstants;

r = alt + Earth.R;

x = r*cos(lat)*cos(lon);
y = r*cos(lat)*sin(lon);
z = r*sin(lat);

r_ECEF = [x;y;z];