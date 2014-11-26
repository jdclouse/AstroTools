function [local_time, ipp_lat, ipp_lon] = ipp_local_time( llh, az, el, UTC )
%ipp_local_time Find the local time of the Ionospheric Pierce Point
%   Detailed explanation goes here
%% Initialize
Re = 6378137.0; % m, USE REAL RE?
ipp_alt = 350e3; % m

%% Find the r of the IPP using vector math
% Vectors in ECEF frame.
lat = llh(1);
lon = llh(2);
r_rec = ellipsoidal2cart(llh(1), llh(2), llh(3)); %?
r_rec_sat = azel2ECEFunitvec(llh, az, el);

% law of sines to find the length of the rec->IPP vector
% NOTE: assuming equatorial radius for all these calculations. The
%       difference between the two is ~0.35% error.
eca = asin(Re*sin(el+pi/2)/(Re+ipp_alt));
tau = pi/2 - el - eca;
l_rec_ipp = Re*sin(tau)/sin(eca);

r_rec_ipp = l_rec_ipp*r_rec_sat/norm(r_rec_sat);
r_ipp = r_rec+r_rec_ipp;

%% Local time using longitude and UTC
[ipp_lat, ipp_lon, ~] = ECEF2ellipsoidal(r_ipp);

local_time = UTC + ipp_lon*43200/pi; %NO

end

