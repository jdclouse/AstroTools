function local_time = ipp_local_time( llh, az, el, UTC )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
%% Initialize
Re = 6378137.0; % m, USE REAL RE?
ipp_alt = 350e3; % m

%% Find the r of the IPP using vector math
% Vectors in ECEF frame.
r_rec = ellipsoidal2cart(llh(1), llh(2), llh(3)); %?
r_rec_sat = azel2ECEFunitvec(llh, az, el);

% law of sines to find the length of the rec->IPP vector
% NOTE: assuming equatorial radius for all these calculations. The
%       difference between the two is ~0.35% error.
int_angle = asin(Re*sin(el+pi/2)/(Re+ipp_alt));
tau = pi/2 - el - int_angle;
l_rec_ipp = Re*sin(tau)/sin(int_angle);

r_rec_ipp = l_rec_ipp*r_rec_sat/norm(r_rec_sat);
r_ipp = r_rec+r_rec_ipp;

%% Local time using longitude and UTC
llh_ipp = ECEF2latlong(r_ipp);
local_time = llh_ipp(2)*UTC; %NO

end

