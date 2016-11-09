function DCM = Euler2DCM( seq_string, angle_vector ) 
%Euler2DCM Turn an Euler Angle set into a DCM
%   Angle vector in radians
% fcnPrintQueue(mfilename('fullpath'))

DCM = eye(3);
%get the trig functions
c = zeros(3,1);
s = zeros(3,1);
c(1) = cos(angle_vector(1));
s(1) = sin(angle_vector(1));
c(2) = cos(angle_vector(2));
s(2) = sin(angle_vector(2));
c(3) = cos(angle_vector(3));
s(3) = sin(angle_vector(3));

for idx = 3:-1:1
    if strcmp(seq_string(idx),'1')
        M = [1 0 0; 0 c(idx) s(idx); 0 -s(idx) c(idx)];
        DCM = DCM*M;
    elseif strcmp(seq_string(idx),'2')
        M = [c(idx) 0 -s(idx); 0 1 0; s(idx) 0 c(idx)];
        DCM = DCM*M;
    elseif strcmp(seq_string(idx),'3')
        M = [c(idx) s(idx) 0; -s(idx) c(idx) 0; 0 0 1];
        DCM = DCM*M;
    else
        fprintf('%s is not a valid axis\n', seq_string(idx))
    end
end
        

end

