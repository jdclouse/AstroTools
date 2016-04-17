function T = gravity_gradient_torque(I,R,params)
% Gravity gradient torque about the principle inertia axes. R is in body
% coords.

R_mag = norm(R);
R5 = R_mag*R_mag*R_mag*R_mag*R_mag;
T = 3*params.G*params.Me/R5 * ...
    [R(2)*R(3)*(I(3)-I(2))...
    R(1)*R(3)*(I(1)-I(3))...
    R(1)*R(2)*(I(2)-I(1))];

