obs = [
6.1773780845922   0                 ;
5.56327661282686  1.31285863495514;
5.69420161397342 -1.54488114381612;
6.15294262127432  0.534923988815733;
5.46251322092491  0.884698415328368;
5.83638064328625 -1.56123248918054;
6.08236452736002  1.00979943157547;
5.40737619817037  0.31705117039215;
5.97065615746125 -1.37453070975606;
5.97369258835895  1.36768169443236;
5.40669060248179 -0.30211158850316];

k1 = 2.5;
k2 = 3.7;
m = 1.5;
h = 5.4;
x0 = 3;
xd0 = 0;
x0 = 4;
xd0 = .2;
last_x_est = Mxx\N;
x0 = 4 + x_est(1);
xd0 = .2 + x_est(2);

% m = 1.6;
% k2 = 3.5;
% h = 5.7;
m = 1.2;
k2 = 3.0;
h = 6.4;

w = sqrt((k1+k2)/m);

STM = @(t) [cos(w*t), sin(w*t)/w;-w*sin(w*t), cos(w*t)];
B1 = x0/2/m;
B2 = xd0/2/m;
get_x = @(t) x0*cos(w*t)+xd0/w*sin(w*t);
get_v = @(t) xd0*cos(w*t)-x0*w*sin(w*t);
Theta = @(t) [B2/w*sin(w*t) + t*(-B2*cos(w*t)+B1*w*sin(w*t)), -(B2/w*sin(w*t) + t*(-B2*cos(w*t)+B1*w*sin(w*t)))/w/w, 0;...
    B1/w*sin(w*t) + t*(B2*w*sin(w*t)+B1*w*w*cos(w*t)), -(B1/w*sin(w*t) + t*(B2*w*sin(w*t)+B1*w*w*cos(w*t)))/w/w, 0];
get_range = @(x) sqrt(x*x+h*h);
get_H_tilde_x = @(x, xdot, rho) [x/rho 0; xdot/rho-x*x*xdot/rho/rho/rho x/rho];
get_H_tilde_c = @(x, xdot, rho) [0 0 h/rho;0 0 -x*xdot*h/rho/rho/rho];

X0_ap = [4;.2] + last_x_est;
Pbar_ap = [100 0;0 10];
c = [.3;.7;-1];
xbar = [x0;xd0]-X0_ap;
Pcc = [0.09 0 0;0 0.49 0; 0 0 1];
% Pcc = [0.01 0 0;0 0.04 0; 0 0 0.09];
R = 1;

Px_bar = Pbar_ap;
Mcc_bar = inv(Pcc);
Mxx_bar = inv(Px_bar);
% Mxc_bar = zeros(1);
Mxc_bar = zeros(2,3);

Mxx = Mxx_bar;
Mxc = Mxc_bar;
Mcc = Mcc_bar;
N = Pbar_ap\xbar;
Hx = [];
Hc = [];
for t = 0:10
    H_tilde_x = get_H_tilde_x(get_x(t), get_v(t), get_range(get_x(t)));
    H_tilde_c = get_H_tilde_c(get_x(t), get_v(t), get_range(get_x(t)));
    Hxi = H_tilde_x*STM(t);
    Hci = H_tilde_x*Theta(t) + H_tilde_c;
    Mxx = Mxx+ (Hxi)'/R*(Hxi);
    Mxc = Mxc+ Hxi'/R*Hci;
    Mcc = Mcc + Hci'/R*Hci;
    
    y = obs(t+1,:)' - [get_range(get_x(t)); get_x(t)*get_v(t)/get_range(get_x(t))]
    
    N = N + Hxi'*R*y;
    
    Hx = [Hx;Hxi];
    Hc = [Hc;Hci];
end

Px = inv(Mxx);
Sxc = -Px*Mxc;

Pxx = Px + Sxc*Pcc*Sxc';
Pxc = Sxc*Pcc;
Gamma = Sxc*diag(sqrt(diag(Pcc)))

Pc = [Pxx Pxc;Pxc' Pcc]
x_est = Mxx\N - Mxx\Mxc*c
