%% Lagrange CR3BP
function state_dot = Lagrange_CR3BP(~, state, opts)

n = opts.n;
M1 = opts.M1;
M2 = opts.M2;
nu = opts.nu;
G = opts.G;
R = opts.R;

x = state(1);
y = state(2);
z = state(3);
x_dot = state(4);
y_dot = state(5);
z_dot = state(6);

r_wrt_1_3 = norm([x+nu*R y z]');
r_wrt_1_3 = r_wrt_1_3*r_wrt_1_3*r_wrt_1_3;
r_wrt_2_3 = norm([x+-(1-nu)*R y z]');
r_wrt_2_3 = r_wrt_2_3*r_wrt_2_3*r_wrt_2_3;

GM1_term = -G*M1/r_wrt_1_3;
GM2_term = -G*M2/r_wrt_2_3;

x_dd = n*n*x + n*y_dot + GM1_term*(x+nu*R) ...
    + GM2_term*(x-(1-nu)*R) +n*y_dot;

y_dd = n*n*y - n*x_dot +GM1_term*y + GM2_term*y -n*x_dot;

z_dd = GM1_term*z + GM2_term*z;

state_dot = [x_dot y_dot z_dot x_dd y_dd z_dd]';