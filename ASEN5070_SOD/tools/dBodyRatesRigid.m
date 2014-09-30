function w_dot = dBodyRatesRigid( w_body, inertia_tensor, cm_torque ) 
%BmatEuler Turn an Euler Angle set into a B matrix
%   theta_dot_vec = B*body_rates_vec
fcnPrintQueue(mfilename('fullpath'))

% right_side = -vecSkew(w_body)*inertia_tensor*w_body + cm_torque;
% vecSkew(w_body)
% -vecSkew(w_body)*inertia_tensor
% inv_I =inv(inertia_tensor);
w_dot = inv(inertia_tensor)*(-vecSkew(w_body)*inertia_tensor*w_body + ...
    cm_torque);
end

