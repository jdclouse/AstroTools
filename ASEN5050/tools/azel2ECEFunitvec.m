function r_hat = azel2ECEFunitvec(llh, az, el)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

% ENU unit vector first
cosel = cos(el);
up = sin(el);
north = cos(az)*cosel;
east = sin(az)*cosel;

latgd = llh(1);
lon = llh(2);
% sines and cosines used later
sinp = sin(latgd);
cosp = cos(latgd);
sinl = sin(lon);
cosl = cos(lon);
% Rotation from ECEF to ENU
R_ecef2local = ...
    [-sinl cosl 0;
    -sinp*cosl -sinp*sinl cosp;
    cosp*cosl cosp*sinl sinp];

r_hat = R_ecef2local'*[east; north; up];

end

