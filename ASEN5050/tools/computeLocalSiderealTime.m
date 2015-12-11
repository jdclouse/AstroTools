function LST = computeLocalSiderealTime(JD, lon, spin_rate)
%computeLocalSiderealTime return local sidereal time
% units in radians
fcnPrintQueue(mfilename('fullpath')) % Add this code to code app 

% Couldn't get 3-45 and 3-46 to work :(
T_UT1 = ((JD)-2451545)/36525;

% GMST0 = 100.4606184 ...
%     + 36000.77005361 * T_UT1 ...
%     + 0.00038793 * T_UT1 * T_UT1 ... 
%     - 2.6e-8 * T_UT1 * T_UT1 * T_UT1;
% GMST0 = GMST0 * pi/180;

GMST = 67310.54841 ...
    + (876600*3600+8640184.812866) * T_UT1 ...
    + 0.093104 * T_UT1 * T_UT1 ...
    - 6.2e-6 * T_UT1 * T_UT1 * T_UT1;
% GMST = GMST * pi/3600/12;
GMST = GMST-86400*floor(GMST/(86400));
GMST = GMST * pi/3600/12;

% GMST = GMST-2*pi*floor(GMST/(2*pi));
% GMST = GMST-360*floor(GMST/(360));

% GMST = GMST0 + spin_rate*(JD - floor(JD));

LST = GMST + lon;