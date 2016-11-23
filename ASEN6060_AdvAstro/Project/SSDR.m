%% SSDR - single shooting differential corrector.
function [X_out, T_out] = SSDR( X_init, T_init, prop_opts )

ode_opts = odeset('RelTol', 1e-8, 'AbsTol', 1e-9);
Xi = X_init;
Ti = T_init;
dx_store = [];
F_store = [];
for ii = 1:20
[T,X] = ode45(@SSDR_deriv,[0 Ti], [Xi; reshape(eye(6),36,1)], ...
    ode_opts, prop_opts);
params_state_dot = Lagrange_CR3BP(0,X(end,1:6)', prop_opts);
F = [Xi - X(end,1:6)'; dot(Xi(4:6),[1;0;1])];

A = [reshape(X(end,7:end),6,6), params_state_dot;...
     0, 0, 0, 1, 0, 1, 0];

dx = A\F;
Xi = Xi + dx(1:6);
Ti = Ti + dx(7);
dx_store = [dx_store, dx];
F_store = [F_store, F];
end
dx_store
F_store
% reshape(X(end,7:end),6,6)
T_out = T(end);
X_out = X(end,1:6)';