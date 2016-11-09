function state_dot = two_body(~, state, opts)

mu = opts.mu;
r = state(1:3);
r_mag = norm(r);
a = state(7:9);

accel = -mu*r/(r_mag*r_mag*r_mag) + a;

state_dot = zeros(9,1);
state_dot(1:3) = state(4:6);
state_dot(4:6) = accel;
