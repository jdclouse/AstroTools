function T = gravity_gradient_torque(I,R_body,params)
% Gravity gradient torque about the principle inertia axes. R is in body
% coords.

R_mag = norm(R_body);
R5 = R_mag*R_mag*R_mag*R_mag*R_mag;
T = 3*params.mu/R5 * ...
    [R_body(2)*R_body(3)*(I(3)-I(2));...
    R_body(1)*R_body(3)*(I(1)-I(3));...
    R_body(1)*R_body(2)*(I(2)-I(1))];

