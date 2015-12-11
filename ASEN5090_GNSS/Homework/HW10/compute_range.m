function [range0, range1, r_gps] = compute_range(eph, PRN, t, userpos)

% fcnPrintQueue(mfilename('fullpath'))
c = 2.99792458e8; %m/s
w = 7292115.0e-11; %rad/s Earth spin rate

GPS_Week = eph(1,19);

% Get the GPS SV vector
[~, r_gps,~,~,~] = broadcast2xv(eph,[GPS_Week t],PRN);
r_gps = r_gps'; % Vertical vectors
R = norm(r_gps-userpos);
range0 = R; % Uncorrected

% Iterate to correct the range for TOF
tol = 1e-6;
oldR = R + 1;
while abs(R-oldR) > tol

Tt = t - R/c;

[~, r_gps_Tt,~,~,~] = broadcast2xv(eph,[GPS_Week Tt],PRN);
r_gps_Tt = r_gps_Tt'; % Vertical vectors

phi = w*(t-Tt);

C = [cos(phi) sin(phi) 0;
    -sin(phi) cos(phi) 0;
    0 0 1];

r_gps = C*r_gps_Tt;

oldR = R;

R = norm(r_gps-userpos);

end

range1 = R;