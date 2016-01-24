function r = ellipsoidal2cart(geod_latitude, geod_longitude, h)
%ellipsoidal2cart Earth ellipsoidal coords to cartesian
% Longitude and Latitude are input radians
% Height h given in km
% fcnPrintQueue(mfilename('fullpath')); % Add this code to code app

%Misra & Enge, Ch 4
a = 6378.1370;%e3; % km
f = 1/298.257223563; %flattening
e = 2*f + f*f;

sin_phi = sin(geod_latitude);
N = a/sqrt(1-e*e*sin_phi*sin_phi);

r = [(N+h)*cos(geod_latitude)*cos(geod_longitude);
    (N+h)*cos(geod_latitude)*sin(geod_longitude);
    (N*(1-e*e)+h)*sin_phi];

end