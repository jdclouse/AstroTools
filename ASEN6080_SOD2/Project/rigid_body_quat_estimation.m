function state_dot = rigid_body_quat_estimation(t, state, opts)
% Assume first 4 states are quaternion
% Next 3 are body rates about principle inertia axes
% Next 3 are principle inertias

q = state(1:4);
w = state(5:7);
I = opts.I;%state(8:10);
R = opts.R;
Rb = EulerParam2DCM(q)*R;
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
B0 = q(1);
B1 = q(2);
B2 = q(3);
B3 = q(4);
dR1_dB = 2*[B0 B3 -B2;
    B1 B2 B3;
    -B2 B1 -B0;
    B3 B0 B1]*R;
dR2_dB = 2*[-B3 B0 B1;
    B2 -B1 B0;
    B1 B2 B3;
    -B0 -B3 B2]*R;
dR3_dB = 2*[B2 -B1 B0;
    B3 -B0 -B1;
    B0 B3 -B2;
    B1 B2 B3]*R;


A = [...
    0 -w1 -w2 -w3 -B1 -B2 -B3;
    w1 0 w3 -w2 B0 -B3 B2;
    w2 -w3 0 w1 B3 B0 -B1;
    w3 w2 -w1 0 -B2 B1 B0;
%     zeros(3,4) eye(3)
    [(R2*dR3_dB' + dR2_dB'*R3)*(I(3)-I(2))/I(1);
    (R1*dR3_dB' + dR1_dB'*R3)*(I(1)-I(3))/I(2);
    (R1*dR2_dB' + dR1_dB'*R2)*(I(2)-I(1))/I(3)]*3*opts.mu/norm(R)^5, ...
    [[0 w3 w2]*-(I(3)-I(2))/I(1);
    [w3 0 w1]*-(I(1)-I(3))/I(2); 
    [w2 w1 0]*-(I(2)-I(1))/I(3)]
    ];

STM = reshape(state(8:end),7,7);
STM_dot = -STM\A;
    
state_dot = [dEP(q,w); w_dot; reshape(STM_dot,49,1)];
    
end