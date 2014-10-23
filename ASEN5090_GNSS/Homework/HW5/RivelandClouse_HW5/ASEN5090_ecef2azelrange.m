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

%Source: Fundamentals of Astrodynamics and Applications 4th Edition
%David Vallado

%Compute range vector
RhoVecEcef = r_sat' - r_site';

%Convert to radians
latgr = latgd * pi/180;
lonr = lon * pi/180;

%Compute Dcm to transform vector from Ecef frame to SEZ frame
DcmEcefToSez = [sin(latgr)*cos(lonr), sin(latgr)*sin(lonr), -cos(latgr); ...
                -sin(lonr), cos(lonr), 0.0; ...
                cos(latgr)*cos(lonr), cos(latgr)*sin(lonr), sin(latgr)];
            
%Compute range vector in the SEZ frame
RhoVecSez = DcmEcefToSez*RhoVecEcef;

%Compute range, magnitude of range vector
range = norm(RhoVecSez);

el = asin(RhoVecSez(3)/range);

if el ~= pi/2
    sinAz = RhoVecSez(2)/(sqrt(RhoVecSez(1)^2 + RhoVecSez(2)^2));
    cosAz = -RhoVecSez(1)/(sqrt(RhoVecSez(1)^2 + RhoVecSez(2)^2));
    
    az = atan2(sinAz, cosAz);
    az = modtwopi(az);
else
    az = 0.0;
end

az = az * 180/pi;
el = el * 180/pi;

