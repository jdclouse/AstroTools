%% SSDR - single shooting differential corrector.
function [X_out, T_out] = SSDR( X_init, T_init, fam_tan, delta_s, prop_opts )

ode_opts = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);
% predict the next orbit w/ pseudo arc len continuation
% delta_s = 1000;
Xi = X_init + delta_s*fam_tan(1:end-1);
Ti = T_init + delta_s*fam_tan(end);
% dx_store = [];
% F_store = [];
% dx_norm_store = [];
% we are starting from a previously-defined orbit
x_tilde_dot = Lagrange_CR3BP(0,X_init, prop_opts);

dx = ones(7,1);

for ii = 1:20
    if norm(dx) < 1e-9
        break;
    end
[T,X] = ode45(@SSDR_deriv,[0 Ti], [Xi; reshape(eye(6),36,1)], ...
    ode_opts, prop_opts);
params_state_dot = Lagrange_CR3BP(0,X(end,1:6)', prop_opts);
% Xi(1)-X_init(1)
curr_delta_s = dot([Xi-X_init;Ti-T_init],fam_tan) ;
% F = [Xi - X(end,1:6)'; dot(Xi-X_init,x_tilde_dot); -delta_s+curr_delta_s];
F = -[X(end,1:6)' - Xi; dot(Xi-X_init,x_tilde_dot); +curr_delta_s-delta_s];

A = [reshape(X(end,7:end),6,6)-eye(6), params_state_dot;...
     x_tilde_dot', 0; fam_tan'];

% dx = A\F;
dx = pinv(A)*F;
Xi = Xi + dx(1:6);
Ti = Ti + dx(7);
% dx_store = [dx_store, dx];
% F_store = [F_store, F];
% dx_norm_store = [dx_norm_store norm(dx)];
end
if ii == 20
    fprintf('Hit 20 loops');
end
% dx_store
% F_store
% dx_norm_store
% reshape(X(end,7:end),6,6)
T_out = T(end);
X_out = X(end,1:6)';