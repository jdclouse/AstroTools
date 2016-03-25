x0 = 1;
w = 10;
l = 1;
g = w*w*l;
theta_0 = 0.01;
theta_dot_0 = 0;
R = 1;

STM = @(t) [cos(w*t), sin(w*t)/w;-w*sin(w*t), cos(w*t)];
Theta = @(t) [0 -t/(2*sqrt(g*l))*(theta_0*sin(w*t) - theta_dot_0/w*cos(w*t));...
    0 -t/(2*sqrt(g*l))*(w*theta_0*cos(w*t) + theta_dot_0*sin(w*t))];
% Theta = @(t) [0 -t/(2*sqrt(g*l))*(theta_0*sin(w*t) - theta_dot_0/w*cos(w*t)) - theta_dot_0/(2*w*g)*sin(w*t);...
%     0 -t/(2*sqrt(g*l))*(w*theta_0*cos(w*t) + theta_dot_0*sin(w*t)) - theta_0/(2*sqrt(g*l))*sin(w*t)];
% 
% Theta = @(t) [0 -theta_0*t/(2*sqrt(l*g))*sin(w*t)+theta_dot_0/w*t/(2*sqrt(l*g))*cos(w*t)-theta_dot_0/2/w/g*sin(w*t);...
%     0 -theta_0/(2*sqrt(l*g))*sin(w*t)-theta_0*w*t/(2*sqrt(l*g))*cos(w*t)-theta_dot_0*t/(2*sqrt(l*g))*sin(w*t)];

get_theta = @(t) theta_0*cos(w*t);
get_range = @(theta) sqrt((x0+l*theta)^2+(l*(1-theta^2/2))^2);
get_H_tilde_x = @(theta, range) [(l^2*theta^3 + 2*x0*l)/2/range, 0];
get_H_tilde_c = @(theta, range) [(2*x0 + 2*l*theta)/2/range, 0];
get_H_tilde_x = @(theta, range) [(l*x0*cos(theta))/range, 0];
get_H_tilde_c = @(theta, range) [(x0 + l*sin(theta))/range, 0];
H_tilde_x = [l 0];
H_tilde_c = [1 0];

% Px_bar  = eye(2);
% Pcc = eye(2);
% Pxc_bar = zeros(2);
% 
% STM = @(t) [1 t; 0 1];
% Theta = @(t) [t*t/2;t];
% 
% H_tilde_x = [l 0];
% H_tilde_c = [0];

Px_bar  = eye(2);
% Pcc = sym('PI');
Pcc = eye(2);
Pxc_bar = zeros(1);

Mcc_bar = inv(Pcc);
Mxx_bar = inv(Px_bar);
% Mxc_bar = zeros(1);
Mxc_bar = zeros(2);

Mxx = Mxx_bar;
Mxc = Mxc_bar;
Mcc = Mcc_bar;
Hx = [];
Hc = [];
for t = 0:2
    H_tilde_x = get_H_tilde_x(get_theta(t),get_range(get_theta(t)));
    H_tilde_c = get_H_tilde_c(get_theta(t),get_range(get_theta(t)));
    Hxi = H_tilde_x*STM(t);
    Hci = H_tilde_x*Theta(t) + H_tilde_c;
    Mxx = Mxx+ (Hxi)'/R*(Hxi);
    Mxc = Mxc+ Hxi'/R*Hci;
    Mcc = Mcc + Hci'/R*Hci;
    
    Hx = [Hx;Hxi];
    Hc = [Hc;Hci];
end

Px = inv(Mxx);
Sxc = -Px*Mxc;

Pxx = Px + Sxc*Pcc*Sxc';
Pxc = Sxc*Pcc;

Pc = [Pxx Pxc;Pxc' Pcc]