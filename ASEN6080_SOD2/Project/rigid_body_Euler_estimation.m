function state_dot = rigid_body_Euler_estimation(t, state, opts)
% Assume first 4 states are quaternion
% Next 3 are body rates about principle inertia axes
% Next 3 are principle inertias

EulerAng = state(1:3);
w = state(4:6);
I = opts.I;%state(8:10);
R = opts.R;
Rb = Euler3212C(EulerAng)*[0;0;norm(R)];
R1 = Rb(1);
R2 = Rb(2);
R3 = Rb(3);
T = zeros(3,1);
if opts.gravity_gradient.use
    T = T + gravity_gradient_torque(I, Rb, opts);
end
w_dot = zeros(3,1);
w_dot(1) = (-(I(3)-I(2))*w(2)*w(3) + T(1))/I(1);
w_dot(2) = (-(I(1)-I(3))*w(1)*w(3) + T(2))/I(2);
w_dot(3) = (-(I(2)-I(1))*w(1)*w(2) + T(3))/I(3);

w1 = w(1);
w2 = w(2);
w3 = w(3);
psi = EulerAng(1);
theta = EulerAng(2);
phi = EulerAng(3);

dR1_dB = [0; -cos(theta); 0]*norm(R);
dR2_dB = [0; -sin(phi)*sin(theta); cos(phi)*cos(theta)]*norm(R);
dR3_dB = [0; -cos(phi)*sin(theta); -sin(phi)*cos(theta)]*norm(R);

Omega = sqrt(opts.mu/norm(R)^3);
w_ON = Omega/cos(theta)*[sin(theta)*sin(psi);-cos(theta)*cos(psi);-sin(psi)];

linearized_correction = ...
    -Omega*[cos(psi)*tan(theta), sec(theta)^2*sin(psi), 0;
    --sin(psi), 0, 0;
    -cos(psi)/cos(theta), -sec(theta)*tan(theta)*sin(psi), 0];
A = [...
    0 sec(theta)*tan(theta)*(sin(phi)*w2 + cos(phi)*w3) sec(theta)*(cos(phi)*w2-sin(phi)*w3) 0 sin(phi)/cos(theta) cos(phi)/cos(theta);
    0 0 -sin(phi)*w2-cos(phi)*w3 0 cos(phi) -sin(phi);
    0 sec(theta)^2*(sin(phi)*w2+cos(phi)*w3) cos(phi)*tan(theta)*w2-sin(phi)*tan(theta)*w3 1 sin(phi)*tan(theta) cos(phi)*tan(theta);
%     zeros(3,4) eye(3)
    [(R2*dR3_dB' + dR2_dB'*R3)*(I(3)-I(2))/I(1);
    (R1*dR3_dB' + dR1_dB'*R3)*(I(1)-I(3))/I(2);
    (R1*dR2_dB' + dR1_dB'*R2)*(I(2)-I(1))/I(3)]*3*opts.mu/norm(R)^5, ...
    [[0 w3 w2]*-(I(3)-I(2))/I(1);
    [w3 0 w1]*-(I(1)-I(3))/I(2); 
    [w2 w1 0]*-(I(2)-I(1))/I(3)]
    ];

A(1:3, 1:3) = A(1:3, 1:3) + linearized_correction;

STM = reshape(state(7:end),6,6);
STM_dot = -STM\A;

state_dot = [dEuler321(EulerAng,w)-w_ON; w_dot; reshape(STM_dot,36,1)];
    
end