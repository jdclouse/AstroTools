function PRV = Euler2PRV( seq_string, angle_vector ) 
%Euler2PRV Turn an Euler Angle set into a PRV
%   PRV = phi*e_vec
fcnPrintQueue(mfilename('fullpath'))

PRV = EulerParam2PRV(Euler2EP(seq_string, angle_vector ));
% DCM = Euler2DCM(seq_string, angle_vector);
% 
% phi = acos((trace(DCM)-1)/2); % gives a phi <= pi
% e = [DCM(2,3)-DCM(3,2); DCM(3,1)-DCM(1,3); DCM(1,2)-DCM(2,1)]/2/sin(phi);
% 
% PRV = phi*e;

end

