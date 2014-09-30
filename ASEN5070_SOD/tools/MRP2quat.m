function quat = MRP2quat( sigma ) 
%BmatEuler Turn an Euler Angle set into a B matrix
%   theta_dot_vec = B*body_rates_vec
fcnPrintQueue(mfilename('fullpath'))

sig_sq=dot(sigma,sigma);

quat = [(1-sig_sq)/(1+sig_sq);...
    2/(1+sig_sq)*sigma];

end

