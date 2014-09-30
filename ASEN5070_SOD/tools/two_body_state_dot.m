function state_dot = two_body_state_dot(t, state, opts)
%two_body_state_dot   Return state_dot given state. Used for numerical
%integration
fcnPrintQueue(mfilename('fullpath')) % Add this code to code app 

state_dot = zeros(6,1);
state_dot(1:3) = state(4:6);
mu = 3.986e5; % km3/s2

r_vec = (state(1:3));
r = norm(r_vec);
state_dot(4:6) = -mu * r_vec/(r*r*r);

if isfield(opts, 'J2')
    if opts.J2 == 1
        state_dot(4:6) = state_dot(4:6) + J2_accel( state(1:3) );
    end
end
if isfield(opts, 'drag')
%     opts.drag.Cd
%     opts.drag.A
%     opts.drag.m
    if opts.drag.use == 1
        state_dot(4:6) = state_dot(4:6) + drag_accel( state, opts.drag );
    end
end