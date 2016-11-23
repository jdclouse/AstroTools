%% Project Find Family
% John Clouse
close all
clear all

% Initialize some things.
addpath('C:\Users\John\Documents\Astro\ASEN5050\tools');
addpath('C:\Users\John\Documents\Astro\ASEN5010\tools');
CelestialConstants;

hw_pub.figWidth = 1120; % pixels
hw_pub.figHeight = 840; % pixels
hw_pub.figPosn = [0, 0, hw_pub.figWidth, hw_pub.figHeight];
% Example: some_fig = figure('Position', hw_pub.figPosn);
hw_pub.lineWidth = 2; % pixels
hw_pub.fontSize = 12;

% ode_opts = odeset('RelTol', 3e-14, 'AbsTol', 1e-20);
ode_opts = odeset('RelTol', 1e-8, 'AbsTol', 1e-9);

% Earth
M1 = Earth.m;

% Moon
M2 = Moon.m;

R = Moon.a;
n = sqrt(G*(M1+M2)/R^3);

moon_dist = R - R*M2/(M1+M2);

propatagor_opts.M1 = M1;
propatagor_opts.M2 = M2;
propatagor_opts.R = R;
propatagor_opts.G = G;
propatagor_opts.nu = M2/(M1+M2);
nu = M2/(M1+M2);
propatagor_opts.n = n;

P = 2*pi/n;

% start with a circular orbit around the moon
alt_init = 1000; %km
r_init = Moon.R+alt_init;
v_init = sqrt(Moon.mu/r_init);
T_init = 2*pi*sqrt(r_init^3/Moon.mu);

% Convert to the synodic frame, starting on the line between Earth and
% Moon.
first_guess_state = [moon_dist-r_init;0;0;0;v_init;0];
% int_time = [0 P/2];
% [Times,X] = ...
%     ode45(@Lagrange_CR3BP,int_time, first_guess_state, ...
%     ode_opts, propatagor_opts);
% 
% figure
% plot(X(:,1),X(:,2))
% hold on
% axis equal
% figure
% plot(X(:,3))
% hold on
% axis equal
d = [1;1];
tol = 1e-13;
X_ini = first_guess_state;
% while abs(d(1)) > tol && abs(d(2)) > tol
while abs(d) > tol 
    X = [X_ini; reshape(eye(6),36,1)];

    [T_out,X_out] = ...
        ode45(@SSDR_deriv, [0,P/2], X, ...
        odeset('Events', @y_crossing),propatagor_opts);

%     d = -[X_out(end,4); X_out(end,6)];
    d = -X_out(end,4);
%     d = [X_ini(1);0]-[X_out(end,1); X_out(end,4)];
    % STM
    STM = reshape(X_out(end,7:end),6,6);
    y_dot = X_out(end,5);
    state_dot = Lagrange_CR3BP(0,X_out(end,1:6)',propatagor_opts);

    % The correction while holding x constant
%     correction = ([STM(4,1) STM(4,5); STM(6,1) STM(6,5)] ...
%         - 1/y_dot*[state_dot(4);state_dot(6)]*[STM(2,1) STM(2,5)])\d;
    correction = d/( STM(4,5) ...
        - 1/y_dot*state_dot(4)* STM(2,5));
%     correction = pinv([STM(1,5); STM(4,5)] ...
%         - 1/y_dot*[state_dot(1);state_dot(4)]*STM(2,5))*d;

%     X_ini(1) = X_ini(1) + correction(1);
%     X_ini(5) = X_ini(5) + correction(2);
    X_ini(5) = X_ini(5) + correction;
    d
end
[Times,X] = ...
    ode45(@Lagrange_CR3BP,[0 10*T_out(end)], X_ini, ...
    ode_opts, propatagor_opts);
figure
plot(X(:,1),X(:,2))
hold on
axis equal


%%

[X_out, T_out] = SSDR(first_guess_state, T_init, propatagor_opts);

[Times,X] = ...
    ode45(@Lagrange_CR3BP,[0 10*T_out], X_out, ...
    ode_opts, propatagor_opts);
figure
plot(X(:,1),X(:,2))
hold on
axis equal