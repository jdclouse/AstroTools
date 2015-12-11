function r_ecef = eci2ecef( r_eci, greenwich_time)
%eci2ecef return ECEF position from ECI coords
% units in radians
fcnPrintQueue(mfilename('fullpath')) % Add this code to code app 

% Just rotate about the z axis.
r_ecef = Euler2DCM('3', [greenwich_time]) * r_eci;