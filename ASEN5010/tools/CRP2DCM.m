function DCM = CRP2DCM( CRP )
%EulerParam2PRV Get PRV from EP

fcnPrintQueue(mfilename('fullpath'))

DCM = ((1-CRP'*CRP)*eye(3) + 2*(CRP*CRP') - 2*vecSkew(CRP))/(1+CRP'*CRP);

end

