function H_tilda = H_tilda_state_rv(state, consts)
%stat_od_proj_H_tilda Calculate H_tilda matrix for Stat OD project
fcnPrintQueue(mfilename('fullpath')) % Add this code to code app 

% Init H_tilda, set up local vars
x = state(1);
y = state(2);
z = state(3);
xdot = state(4);
ydot = state(5);
zdot = state(6);

theta_dot = consts.theta_dot;
theta = consts.t*consts.theta_dot;

xs = consts.site_r(1);
ys = consts.site_r(2);
zs = consts.site_r(3);

H_tilda = zeros(2,6);
H_tilda(1,:) = [                                                                                                                                                                                                                                                                                                                                                    (2*x - 2*xs*cos(theta) + 2*ys*sin(theta))/(2*((ys*cos(theta) - y + xs*sin(theta))^2 + (z - zs)^2 + (x - xs*cos(theta) + ys*sin(theta))^2)^(1/2)),                                                                                                                                                                                                                                                                                                                                                   -(2*ys*cos(theta) - 2*y + 2*xs*sin(theta))/(2*((ys*cos(theta) - y + xs*sin(theta))^2 + (z - zs)^2 + (x - xs*cos(theta) + ys*sin(theta))^2)^(1/2)),                                                                                                                                                                                                                                                                                              (2*z - 2*zs)/(2*((ys*cos(theta) - y + xs*sin(theta))^2 + (z - zs)^2 + (x - xs*cos(theta) + ys*sin(theta))^2)^(1/2)),                                                                                                                                      0,                                                                                                                                       0,                                                                                                           0];
H_tilda(2,:) = [ (xdot + theta_dot*ys*cos(theta) + theta_dot*xs*sin(theta))/((ys*cos(theta) - y + xs*sin(theta))^2 + (z - zs)^2 + (x - xs*cos(theta) + ys*sin(theta))^2)^(1/2) - ((2*x - 2*xs*cos(theta) + 2*ys*sin(theta))*(sin(theta)*(xdot*ys - xs*ydot) - cos(theta)*(xdot*xs + ydot*ys) + x*xdot + y*ydot + z*zdot - zdot*zs + theta_dot*cos(theta)*(x*ys - xs*y) + theta_dot*sin(theta)*(x*xs + y*ys)))/(2*((ys*cos(theta) - y + xs*sin(theta))^2 + (z - zs)^2 + (x - xs*cos(theta) + ys*sin(theta))^2)^(3/2)), (ydot - theta_dot*xs*cos(theta) + theta_dot*ys*sin(theta))/((ys*cos(theta) - y + xs*sin(theta))^2 + (z - zs)^2 + (x - xs*cos(theta) + ys*sin(theta))^2)^(1/2) + ((2*ys*cos(theta) - 2*y + 2*xs*sin(theta))*(sin(theta)*(xdot*ys - xs*ydot) - cos(theta)*(xdot*xs + ydot*ys) + x*xdot + y*ydot + z*zdot - zdot*zs + theta_dot*cos(theta)*(x*ys - xs*y) + theta_dot*sin(theta)*(x*xs + y*ys)))/(2*((ys*cos(theta) - y + xs*sin(theta))^2 + (z - zs)^2 + (x - xs*cos(theta) + ys*sin(theta))^2)^(3/2)), zdot/((ys*cos(theta) - y + xs*sin(theta))^2 + (z - zs)^2 + (x - xs*cos(theta) + ys*sin(theta))^2)^(1/2) - ((2*z - 2*zs)*(sin(theta)*(xdot*ys - xs*ydot) - cos(theta)*(xdot*xs + ydot*ys) + x*xdot + y*ydot + z*zdot - zdot*zs + theta_dot*cos(theta)*(x*ys - xs*y) + theta_dot*sin(theta)*(x*xs + y*ys)))/(2*((ys*cos(theta) - y + xs*sin(theta))^2 + (z - zs)^2 + (x - xs*cos(theta) + ys*sin(theta))^2)^(3/2)), (x - xs*cos(theta) + ys*sin(theta))/((ys*cos(theta) - y + xs*sin(theta))^2 + (z - zs)^2 + (x - xs*cos(theta) + ys*sin(theta))^2)^(1/2), -(ys*cos(theta) - y + xs*sin(theta))/((ys*cos(theta) - y + xs*sin(theta))^2 + (z - zs)^2 + (x - xs*cos(theta) + ys*sin(theta))^2)^(1/2), (z - zs)/((ys*cos(theta) - y + xs*sin(theta))^2 + (z - zs)^2 + (x - xs*cos(theta) + ys*sin(theta))^2)^(1/2)];
