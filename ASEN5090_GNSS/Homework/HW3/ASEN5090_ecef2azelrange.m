function [az,el,range] = ASEN5090_ecef2azelrange(r_sat,r_site,latgd,lon)

%==========================================================================
%==========================================================================
% [az,el,range] = ecef2azelrange(r_sat,r_site,latgd,lon)
%
% Calculates the azimuth, elevation, and range of a satellite with respect
%  to an observation site.
%
%
% Author: Ben K. Bradley
% Date: 11/15/2010
% Modified to remove calculations for ASEN5090 assignments
%
% INPUT:         Description                                         Units
%
%  r_sat      - position of satellite in ECEF frame                 [x y z]
%  r_site     - position of observing site in ECEF frame            [x y z]
%  latgd      - geodetic latitude of observation site          [-90,90] deg
%  lon        - longitude of observation site     [-180,180] or [0,360] deg
%
%
% OUTPUT:       
%    
%  az         - azimuth (degrees clockwise from North)          [0,360] deg
%  el         - elevation (degrees up from horizon)            [-90,90] deg
%  range      - distance from observation site to satellite                                    
%
%
% Coupling:
%
%  none
%
%
%==========================================================================
%==========================================================================

% Satellite pos rel to site
r_site2sat_ecef = r_sat - r_site;
% sines and cosines used later
sinp = sind(latgd);
cosp = cosd(latgd);
sinl = sind(lon);
cosl = cosd(lon);
% Rotation from ECEF to ENU
R_ecef2local = ...
    [-sinl cosl 0;
    -sinp*cosl -sinp*sinl cosp;
    cosp*cosl cosp*sinl sinp];

% Rotate the rel pos into ENU
r_site2sat_enu =  R_ecef2local*r_site2sat_ecef;

% ENU coords give you az/el/range
az = atan2(r_site2sat_enu(1), r_site2sat_enu(2))*180/pi;
range = norm(r_site2sat_enu);
el = asin(r_site2sat_enu(3)/range)*180/pi;
end

