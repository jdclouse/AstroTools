function r_eci = ecef2eci( r_ecef, greenwich_time)
%ecef2eci return ECI position from ECEF coords
% units in radians
fcnPrintQueue(mfilename('fullpath')) % Add this code to code app 

% Just a negative rotation about the z axis.
r_eci = Euler2DCM('3', [-greenwich_time]) * r_ecef;