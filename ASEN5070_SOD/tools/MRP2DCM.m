function DCM = MRP2DCM( sigma ) 
%BmatEuler Turn an Euler Angle set into a B matrix
%   theta_dot_vec = B*body_rates_vec
fcnPrintQueue(mfilename('fullpath'))
skew = vecSkew(sigma);
sig_sq=dot(sigma,sigma);
d = (1+sig_sq)*(1+sig_sq);
DCM = eye(3)+(8*skew*skew - 4*(1-sig_sq)*skew)/d;
end

