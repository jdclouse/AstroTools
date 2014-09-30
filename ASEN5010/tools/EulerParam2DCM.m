function DCM = EulerParam2DCM( EP )
%EulerParam2PRV Get PRV from EP

fcnPrintQueue(mfilename('fullpath'))

B0 = EP(1);
B1 = EP(2);
B2 = EP(3);
B3 = EP(4);

C11 = B0*B0 + B1*B1 - B2*B2 - B3*B3;
C22 = B0*B0 - B1*B1 + B2*B2 - B3*B3;
C33 = B0*B0 - B1*B1 - B2*B2 + B3*B3;
C12 = 2*(B1*B2 + B0*B3);
C13 = 2*(B1*B3 - B0*B2);
C21 = 2*(B1*B2 - B0*B3);
C23 = 2*(B2*B3 + B0*B1);
C31 = 2*(B1*B3 + B0*B2);
C32 = 2*(B2*B3 - B0*B1);
DCM = [C11 C12 C13; C21 C22 C23; C31 C32 C33];

end

