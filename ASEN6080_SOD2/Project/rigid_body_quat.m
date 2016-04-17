function state_dot = rigid_body_quat(state, t)
% Assume first 4 states are quaternion
% Next 3 are body rates about principle inertia axes
% Next 3 are principle inertias

w = state(4:7);
I = state(8:10);
T = zeros(3,1);
w_dot = zeros(3,1);
w_dot(1) = (-(I(3)-I(2))*w(2)*w(3) + T(1))/I(1);
w_dot(2) = (-(I(1)-I(3))*w(1)*w(3) + T(2))/I(2);
w_dot(3) = (-(I(2)-I(1))*w(1)*w(2) + T(3))/I(3);

state_dot = [q_dot w_dot zeros(3,1)];
end