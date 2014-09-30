function DCM = quat2DCM( quat )
fcnPrintQueue(mfilename('fullpath'))
% Scalar term first
B0 = quat(1);
B1 = quat(2);
B2 = quat(3);
B3 = quat(4);

DCM = zeros(3);
DCM(1,1) = B0 + B1 - B2 - B3;
DCM(2,2) = B0 - B1 + B2 - B3;
DCM(3,3) = B0 - B1 - B2 + B3;
DCM(1,2) =  2*(B1*B2 + B0*B3);
DCM(2,1) =  2*(B1*B2 - B0*B3);
DCM(1,3) =  2*(B1*B3 - B0*B2);
DCM(3,1) =  2*(B1*B3 + B0*B2);
DCM(2,3) =  2*(B2*B3 + B0*B1);
DCM(3,2) =  2*(B2*B3 - B0*B1);


end

