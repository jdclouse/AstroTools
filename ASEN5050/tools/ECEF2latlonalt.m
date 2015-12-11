function [lat, lon, alt] = ECEF2latlonalt( r_ecef)
%ECEF2latlonalt return geocentric latitude, longitude, and altitude 
% from ECEF coords.
% units in radians
fcnPrintQueue(mfilename('fullpath')) % Add this code to code app 

% Get some useful constants
CelestialConstants;

% Just subtract the earth radius to get altitude
r = norm(r_ecef);
alt = r - Earth.R;

lat = asin(r_ecef(3)/r);
% y/x = tan(lon)
lon = atan2(r_ecef(2), r_ecef(1));