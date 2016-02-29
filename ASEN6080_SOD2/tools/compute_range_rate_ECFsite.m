function range_rate = compute_range_rate_ECFsite( inrtl_state, ...
    ecf_site, theta, theta_dot )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
% fcnPrintQueue(mfilename('fullpath'))

x = inrtl_state(1);
y = inrtl_state(2);
z = inrtl_state(3);
xdot = inrtl_state(4);
ydot = inrtl_state(5);
zdot = inrtl_state(6);
xs = ecf_site(1);
ys = ecf_site(2);
zs = ecf_site(3);

range_rate = (x*xdot+y*ydot+z*zdot...
    - (xdot*xs+ydot*ys)*cos(theta) + theta_dot*(x*xs + y*ys)*sin(theta) ...
    + (xdot*ys-ydot*xs)*sin(theta) + theta_dot*(x*ys-y*xs)*cos(theta) ...
    - zdot*zs)/compute_range_ECFsite(inrtl_state(1:3), ecf_site, theta);

end

