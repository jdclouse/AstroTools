function euler_param_dot = dEulerParam( EP, w_body ) 
%BmatEuler Turn an Euler Angle set into a B matrix
%   theta_dot_vec = B*body_rates_vec
fcnPrintQueue(mfilename('fullpath'))
B0=EP(1);
B1=EP(2);
B2=EP(3);
B3=EP(4);
euler_param_dot = 0.5*...
    [B0, -B1, -B2, -B3;...
    B1, B0, -B3, B2;...
    B2, B3, B0, -B1;...
    B3, -B2, B1, B0] * [0;w_body];
end

