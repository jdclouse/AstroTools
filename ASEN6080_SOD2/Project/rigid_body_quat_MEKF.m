function state_dot = rigid_body_quat_MEKF(t, state, opts)
% Assume first 4 states are quaternion
% Next 3 are body rates about principle inertia axes
% Next 3 are principle inertias

q = state(1:4);
w = state(5:7);
I = opts.I;%state(8:10);
% R = opts.R;
% Rb = EulerParam2DCM(q)*R;
T = zeros(3,1);
if opts.gravity_gradient.use
    T = T + gravity_gradient_torque(I, Rb, opts);
end
w_dot = zeros(3,1);
% w_dot(1) = (-(I(3)-I(2))*w(2)*w(3) + T(1))/I(1);
% w_dot(2) = (-(I(1)-I(3))*w(1)*w(3) + T(2))/I(2);
% w_dot(3) = (-(I(2)-I(1))*w(1)*w(2) + T(3))/I(3);

% F = [-vec_skew(w) -eye(3); zeros(3,6)];
% G = [-eye(3) zeros(3);zeros(3) eye(3)];
F = -vecSkew(w);
G = eye(3);

P = reshape(state(8:end),3,3);
R = opts.R;
H = eye(3);
Q = opts.Q;

P_dot = F*P + P*F' + G*Q*G -P*H'/R*H*P;
    
state_dot = [dEP(q,w); w_dot; reshape(P_dot,9,1)];
    
end