function B_inv = BinvEuler( seq_string, angle_vector ) 
%BinvEuler Turn an Euler Angle set into a inv(B) matrix
%   theta_dot_vec = B*body_rates_vec
fcnPrintQueue(mfilename('fullpath'))
B_inv = zeros(3);
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
    B_inv = [-s(2), 0, 1; c(2)*s(3), c(3), 0; c(2)*c(3), - s(3), 0];
else
    fprintf('this rotation sequence is not supported');
end
end

