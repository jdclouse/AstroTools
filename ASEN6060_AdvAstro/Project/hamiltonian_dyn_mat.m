function dyn_mat = hamiltonian_dyn_mat(state, opts)

n = opts.n;
M1 = opts.M1;
M2 = opts.M2;
nu = opts.nu;
G = opts.G;
R = opts.R;

x = state(1);
y = state(2);
z = state(3);
% x_dot = state(4);
% y_dot = state(5);
% z_dot = state(6);

r_wrt_1_3 = norm([x+nu*R y z]')^3;
r_wrt_2_3 = norm([x+-(1-nu)*R y z]')^3;
r_wrt_1_5 = norm([x+nu*R y z]')^5;
r_wrt_2_5 = norm([x+-(1-nu)*R y z]')^5;

GM1_term_3 = G*M1/r_wrt_1_3;
GM2_term_3 = G*M2/r_wrt_2_3;
GM1_term_5 = -3*G*M1/r_wrt_1_5;
GM2_term_5 = -3*G*M2/r_wrt_2_5;

dx_dx_H = GM1_term_5*(x+nu*R)*(x+nu*R) ...
    + GM2_term_5*(x-(1-nu)*R)*(x-(1-nu)*R)...
    + GM1_term_3 + GM2_term_3;
dy_dx_H = GM1_term_5*(x+nu*R)*y ...
    + GM2_term_5*(x-(1-nu)*R)*y;
dz_dx_H = GM1_term_5*(x+nu*R)*z ...
    + GM2_term_5*(x-(1-nu)*R)*z;
dp_dx_H = [0 -n 0];
dy_dy_H = GM1_term_5*y*y ...
    + GM2_term_5*y*y...
    + GM1_term_3 + GM2_term_3;
dz_dy_H = GM1_term_5*y*z ...
    + GM2_term_5*y*z;
dp_dy_H = [-n 0 0];
dz_dz_H = GM1_term_5*z*z ...
    + GM2_term_5*z*z...
    + GM1_term_3 + GM2_term_3;
dp_dz_H = [0 0 0];
dp_dp_H = eye(3);

dyn_mat = [zeros(3) eye(3);-eye(3) zeros(3)]*...
    [dx_dx_H dy_dx_H dz_dx_H dp_dx_H;...
     dy_dx_H dy_dy_H dz_dy_H dp_dy_H;...
     dz_dx_H dz_dy_H dz_dz_H dp_dz_H;...
     dp_dx_H' dp_dy_H' dp_dz_H' dp_dp_H];