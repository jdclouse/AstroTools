function state_dot = SSDR_deriv(~, state, opts)
% Uses the same parameters as the CR3BP method
n = opts.n;
M1 = opts.M1;
M2 = opts.M2;
nu = opts.nu;
G = opts.G;
R = opts.R;

STM = reshape(state(7:end),6,6);

x = state(1);
y = state(2);
z = state(3);
x_dot = state(4);
y_dot = state(5);
z_dot = state(6);

params_state_dot = Lagrange_CR3BP(0,[x;y;z;x_dot;y_dot;z_dot], opts);

r_wrt_1_3 = norm([x+nu*R y z]')^3;
r_wrt_2_3 = norm([x+-(1-nu)*R y z]')^3;
r_wrt_1_5 = norm([x+nu*R y z]')^5;
r_wrt_2_5 = norm([x+-(1-nu)*R y z]')^5;

GM1_term_3 = -G*M1/r_wrt_1_3;
GM2_term_3 = -G*M2/r_wrt_2_3;
GM1_term_5 = 3*G*M1/r_wrt_1_5;
GM2_term_5 = 3*G*M2/r_wrt_2_5;

d_xdd_dx = n*n + GM1_term_3 + GM2_term_3 ...
    + GM1_term_5*(x+nu*R)*(x+nu*R)...
    + GM2_term_5*(x-(1-nu)*R)*(x-(1-nu)*R);
d_xdd_dy = GM1_term_5*(x+nu*R)*y + GM2_term_5*(x-(1-nu)*R)*y;
d_xdd_dz = GM1_term_5*(x+nu*R)*z + GM2_term_5*(x-(1-nu)*R)*z;

d_xdd_dqd = [0, 2*n, 0];

d_ydd_dx = GM1_term_5*(x+nu*R)*y + GM2_term_5*(x-(1-nu)*R)*y;
d_ydd_dy = n*n + GM1_term_3 + GM2_term_3 + GM1_term_5*y*y + GM2_term_5*y*y;
d_ydd_dz = GM1_term_5*y*z + GM2_term_5*y*z;

d_ydd_dqd = [-2*n, 0, 0];

d_zdd_dx = GM1_term_5*(x+nu*R)*z + GM2_term_5*(x-(1-nu)*R)*z;
d_zdd_dy = GM1_term_5*y*z + GM2_term_5*y*z;
d_zdd_dz = GM1_term_3 + GM2_term_3 + GM1_term_5*z*z + GM2_term_5*z*z;

A = [zeros(3,3), eye(3);...
    d_xdd_dx, d_xdd_dy, d_xdd_dz, d_xdd_dqd;...
    d_ydd_dx, d_ydd_dy, d_ydd_dz, d_ydd_dqd;...
    d_zdd_dx, d_zdd_dy, d_zdd_dz, 0,0,0];
    
STM_dot = A*STM;
state_dot = [params_state_dot;reshape(STM_dot,36,1)];