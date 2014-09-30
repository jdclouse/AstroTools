function DCM = PRV2DCM( PRV ) 
%PRV2DCM Turn a PRV into a DCM
%   PRV = phi*e
fcnPrintQueue(mfilename('fullpath'))

DCM = zeros(3);
phi = norm(PRV);
e = PRV/phi;
cphi = cos(phi);
sphi = sin(phi);
sigma = 1-cphi;

DCM(1,1) = e(1)*e(1)*sigma + cphi;
DCM(1,2) = e(1)*e(2)*sigma + e(3)*sphi;
DCM(1,3) = e(1)*e(3)*sigma - e(2)*sphi;
DCM(2,1) = e(2)*e(1)*sigma - e(3)*sphi;
DCM(2,2) = e(2)*e(2)*sigma + cphi;
DCM(2,3) = e(2)*e(3)*sigma + e(1)*sphi;
DCM(3,1) = e(3)*e(1)*sigma + e(2)*sphi;
DCM(3,2) = e(3)*e(2)*sigma - e(1)*sphi;
DCM(3,3) = e(3)*e(3)*sigma + cphi;

end

