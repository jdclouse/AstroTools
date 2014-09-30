function PRV = EulerParam2PRV( EP )
%EulerParam2PRV Get PRV from EP

fcnPrintQueue(mfilename('fullpath'))

phi = acos(EP(1))*2;
e = [EP(2); EP(3); EP(4)]/sin(phi/2);

PRV = phi*e;

end

