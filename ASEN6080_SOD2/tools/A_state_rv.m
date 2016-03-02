function A = A_state_rv(state, consts)
%stat_od_proj_A Calculate A matrix for Stat OD project
% fcnPrintQueue(mfilename('fullpath')) % Add this code to code app 

% Init A, set up local vars
% A = zeros(consts.state_len);
x = state(1);
y = state(2);
z = state(3);
xdot = state(4);
ydot = state(5);
zdot = state(6);
mu = consts.mu;
J2 = consts.J2.params.J2;
Cd = 0;

Re = consts.Re;
area = consts.area;
rho = consts.rho;
theta_dot = consts.theta_dot;
m = consts.m;

% vars to reduce computations
x2 = x*x;
y2 = y*y;
z2 = z*z;
r = sqrt(x2+y2+z2);
sqrt_r = sqrt(r);
v = sqrt(xdot*xdot+ydot*ydot+zdot*zdot);
rel_wind_x = (xdot + theta_dot*y);
rel_wind_y = (ydot - theta_dot*x);
zdot2 = zdot*zdot;
rel_wind_mag = sqrt(rel_wind_x*rel_wind_x + rel_wind_y*rel_wind_y + zdot2);
Re2 = Re*Re;

rho0 = 3.614e-13; %kg/m3
r0 = 700000+6378136.3; %km
H = 88667.0; %km
r_sqrd = x^2 + y^2 + z^2;
r_sqrd_3_2 = r_sqrd^(3/2);
r_sqrd_5_2 = r_sqrd^(5/2);
r_sqrd_7_2 = r_sqrd^(7/2);
r_sqrd_9_2 = r_sqrd^(9/2);

% Only a few elements are populated
% A1 = [ 0,0,0,1,0,0];
% A2 = [ 0,0,0,0,1,0];
% A3 = [ 0,0,0,0,0,1];
% A4 = [            (3*mu*x^2)/(x^2 + y^2 + z^2)^(5/2) - mu/(x^2 + y^2 + z^2)^(3/2) + (3*J2*Re^2*mu*((5*z^2)/(x^2 + y^2 + z^2) - 1))/(2*(x^2 + y^2 + z^2)^(5/2)) - (15*J2*Re^2*mu*x^2*z^2)/(x^2 + y^2 + z^2)^(9/2) - (15*J2*Re^2*mu*x^2*((5*z^2)/(x^2 + y^2 + z^2) - 1))/(2*(x^2 + y^2 + z^2)^(7/2)) + (0)/(1) + (0)/(1), (3*mu*x*y)/(x^2 + y^2 + z^2)^(5/2) - (0)/(1) - (15*J2*Re^2*mu*x*y*z^2)/(x^2 + y^2 + z^2)^(9/2) - (15*J2*Re^2*mu*x*y*((5*z^2)/(x^2 + y^2 + z^2) - 1))/(2*(x^2 + y^2 + z^2)^(7/2)) - (0)/(1) + (0)/(1),                                                                                           (3*mu*x*z)/(x^2 + y^2 + z^2)^(5/2) + (3*J2*Re^2*mu*x*((10*z)/(x^2 + y^2 + z^2) - (10*z^3)/(x^2 + y^2 + z^2)^2))/(2*(x^2 + y^2 + z^2)^(5/2)) - (15*J2*Re^2*mu*x*z*((5*z^2)/(x^2 + y^2 + z^2) - 1))/(2*(x^2 + y^2 + z^2)^(7/2)) + (0)/(1), - (0)/(1) - (0)/(1),                                                                                                                                -(0)/(1),                                                                                                             -(0)/(1)];
% A5 = [ (3*mu*x*y)/(x^2 + y^2 + z^2)^(5/2) + (0)/(1) - (15*J2*Re^2*mu*x*y*z^2)/(x^2 + y^2 + z^2)^(9/2) - (15*J2*Re^2*mu*x*y*((5*z^2)/(x^2 + y^2 + z^2) - 1))/(2*(x^2 + y^2 + z^2)^(7/2)) + (0)/(1) + (0)/(1),            (3*mu*y^2)/(x^2 + y^2 + z^2)^(5/2) - mu/(x^2 + y^2 + z^2)^(3/2) + (3*J2*Re^2*mu*((5*z^2)/(x^2 + y^2 + z^2) - 1))/(2*(x^2 + y^2 + z^2)^(5/2)) - (15*J2*Re^2*mu*y^2*z^2)/(x^2 + y^2 + z^2)^(9/2) - (15*J2*Re^2*mu*y^2*((5*z^2)/(x^2 + y^2 + z^2) - 1))/(2*(x^2 + y^2 + z^2)^(7/2)) - (0)/(1) + (0)/(1),                                                                                           (3*mu*y*z)/(x^2 + y^2 + z^2)^(5/2) + (3*J2*Re^2*mu*y*((10*z)/(x^2 + y^2 + z^2) - (10*z^3)/(x^2 + y^2 + z^2)^2))/(2*(x^2 + y^2 + z^2)^(5/2)) - (15*J2*Re^2*mu*y*z*((5*z^2)/(x^2 + y^2 + z^2) - 1))/(2*(x^2 + y^2 + z^2)^(7/2)) + (0)/(1),                                                                                                                                -(0)/(1), - (0)/(1) - (0)/(1),                                                                                                             -(0)/(1)];
% A6 = [                                                                                                                                                        (3*mu*x*z)/(x^2 + y^2 + z^2)^(5/2) - (15*J2*Re^2*mu*x*z^3)/(x^2 + y^2 + z^2)^(9/2) - (15*J2*Re^2*mu*x*z*((5*z^2)/(x^2 + y^2 + z^2) - 3))/(2*(x^2 + y^2 + z^2)^(7/2)) + (0)/(1) + (0)/(1),                                                                                                                                                        (3*mu*y*z)/(x^2 + y^2 + z^2)^(5/2) - (15*J2*Re^2*mu*y*z^3)/(x^2 + y^2 + z^2)^(9/2) - (15*J2*Re^2*mu*y*z*((5*z^2)/(x^2 + y^2 + z^2) - 3))/(2*(x^2 + y^2 + z^2)^(7/2)) - (0)/(1) + (0)/(1), (3*mu*z^2)/(x^2 + y^2 + z^2)^(5/2) - mu/(x^2 + y^2 + z^2)^(3/2) + (3*J2*Re^2*mu*((5*z^2)/(x^2 + y^2 + z^2) - 3))/(2*(x^2 + y^2 + z^2)^(5/2)) + (3*J2*Re^2*mu*z*((10*z)/(x^2 + y^2 + z^2) - (10*z^3)/(x^2 + y^2 + z^2)^2))/(2*(x^2 + y^2 + z^2)^(5/2)) - (15*J2*Re^2*mu*z^2*((5*z^2)/(x^2 + y^2 + z^2) - 3))/(2*(x^2 + y^2 + z^2)^(7/2)) + (0)/(1),                                                                                                                                                -(0)/(1),                                                                                                                                                -(0)/(1), - (0)/(1) - (0)/(1)];
A1 = [ 0,0,0,1,0,0];
A2 = [ 0,0,0,0,1,0];
A3 = [ 0,0,0,0,0,1];
A4 = [            (3*mu*x^2)/(r_sqrd_5_2) - mu/(r_sqrd_3_2) + (3*J2*Re^2*mu*((5*z^2)/(r_sqrd) - 1))/(2*(r_sqrd_5_2)) - (15*J2*Re^2*mu*x^2*z^2)/(r_sqrd_9_2) - (15*J2*Re^2*mu*x^2*((5*z^2)/(r_sqrd) - 1))/(2*(r_sqrd_7_2)) + 0 + 0, (3*mu*x*y)/(r_sqrd_5_2) - 0 - (15*J2*Re^2*mu*x*y*z^2)/(r_sqrd_9_2) - (15*J2*Re^2*mu*x*y*((5*z^2)/(r_sqrd) - 1))/(2*(r_sqrd_7_2)) - 0 + 0,                                                                                           (3*mu*x*z)/(r_sqrd_5_2) + (3*J2*Re^2*mu*x*((10*z)/(r_sqrd) - (10*z^3)/(r_sqrd)^2))/(2*(r_sqrd_5_2)) - (15*J2*Re^2*mu*x*z*((5*z^2)/(r_sqrd) - 1))/(2*(r_sqrd_7_2)), 0, 0, 0];
A5 = [ (3*mu*x*y)/(r_sqrd_5_2) + 0 - (15*J2*Re^2*mu*x*y*z^2)/(r_sqrd_9_2) - (15*J2*Re^2*mu*x*y*((5*z^2)/(r_sqrd) - 1))/(2*(r_sqrd_7_2)) + 0 + 0,            (3*mu*y^2)/(r_sqrd_5_2) - mu/(r_sqrd_3_2) + (3*J2*Re^2*mu*((5*z^2)/(r_sqrd) - 1))/(2*(r_sqrd_5_2)) - (15*J2*Re^2*mu*y^2*z^2)/(r_sqrd_9_2) - (15*J2*Re^2*mu*y^2*((5*z^2)/(r_sqrd) - 1))/(2*(r_sqrd_7_2)) - 0 + 0,                                                                                           (3*mu*y*z)/(r_sqrd_5_2) + (3*J2*Re^2*mu*y*((10*z)/(r_sqrd) - (10*z^3)/(r_sqrd)^2))/(2*(r_sqrd_5_2)) - (15*J2*Re^2*mu*y*z*((5*z^2)/(r_sqrd) - 1))/(2*(r_sqrd_7_2)), 0, 0, 0];
A6 = [ (3*mu*x*z)/(r_sqrd_5_2) - (15*J2*Re^2*mu*x*z^3)/(r_sqrd_9_2) - (15*J2*Re^2*mu*x*z*((5*z^2)/(r_sqrd) - 3))/(2*(r_sqrd_7_2)) + 0 + 0,                                                                                                                                                        (3*mu*y*z)/(r_sqrd_5_2) - (15*J2*Re^2*mu*y*z^3)/(r_sqrd_9_2) - (15*J2*Re^2*mu*y*z*((5*z^2)/(r_sqrd) - 3))/(2*(r_sqrd_7_2)) - 0 + 0, (3*mu*z^2)/(r_sqrd_5_2) - mu/(r_sqrd_3_2) + (3*J2*Re^2*mu*((5*z^2)/(r_sqrd) - 3))/(2*(r_sqrd_5_2)) + (3*J2*Re^2*mu*z*((10*z)/(r_sqrd) - (10*z^3)/(r_sqrd)^2))/(2*(r_sqrd_5_2)) - (15*J2*Re^2*mu*z^2*((5*z^2)/(r_sqrd) - 3))/(2*(r_sqrd_7_2)), -0, -0, 0];



A(1,:) = A1;
A(2,:) = A2;
A(3,:) = A3;
A(4,:) = A4;
A(5,:) = A5;
A(6,:) = A6;
A = [ A1; A2; A3; A4; A5; A6];