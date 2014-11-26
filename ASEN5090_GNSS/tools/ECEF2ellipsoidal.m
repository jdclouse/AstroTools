function [glat, lon, h] = ECEF2ellipsoidal(r)
%ECEF2ellipsoidal Earth-Centered Earth-Fixed coords to ellipsoidal
% Longitude and Latitude are output radians
% Height h given in m


% fcnPrintQueue(mfilename('fullpath')) % Add this code to code app

%Misra & Enge, Ch 4
a = 6378137.0; % km
f = 1/298.257223563; %flattening
e = 2*f + f*f;

x = r(1);
y = r(2);
z = r(3);
lon = atan2(y,x);

p = sqrt(x*x+y*y);

glat = atan2(p,z);

lattol = 0.0001;
htol = 0.0001;
lat_err = lattol+1;
h_err = htol+1;
h = 0;
counter = 0;
while lat_err > lattol && h_err > htol && counter < 5
    hlast = h;
    latlast = glat;
    N = a/sqrt(1-e*e*sin(glat)*sin(glat));
    h = p/cos(glat)-N;
    glat = atan(z/(p*(1-e*e*N/(N+h))));
    lat_err = abs(latlast-glat);
    h_err = abs(hlast-h);
    counter = counter + 1;
end
% 
% llh = [glat lon h]';

end