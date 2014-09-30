function euler_dot = dEuler( seq_string, angle_vector, w_body ) 
%BmatEuler Turn an Euler Angle set into a B matrix
%   theta_dot_vec = B*body_rates_vec
fcnPrintQueue(mfilename('fullpath'))

euler_dot = BmatEuler(seq_string, angle_vector)*w_body;
end

