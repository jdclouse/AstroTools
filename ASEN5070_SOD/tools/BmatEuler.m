function B = BmatEuler( seq_string, angle_vector ) 
%BmatEuler Turn an Euler Angle set into a B matrix
%   theta_dot_vec = B*body_rates_vec
fcnPrintQueue(mfilename('fullpath'))
B = zeros(3);
%get the trig functions
c = zeros(3,1);
s = zeros(3,1);
c(1) = cos(angle_vector(1));
s(1) = sin(angle_vector(1));
c(2) = cos(angle_vector(2));
s(2) = sin(angle_vector(2));
c(3) = cos(angle_vector(3));
s(3) = sin(angle_vector(3));

if strcmp(seq_string,'321')
    B = [0, s(3), c(3); 0, c(2)*c(3), -c(2)*s(3); c(2), s(2)*s(3), s(2)*c(3)] / c(2);
else
    fprintf('this rotation sequence is not supported');
end
end

