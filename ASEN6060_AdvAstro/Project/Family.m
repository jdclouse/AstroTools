%% Project Find Family
% John Clouse
close all
clear all

% Initialize some things.
% root = 'C:\Users\John\Documents\Astro';
root = 'C:\Users\jclouse1\Documents\Class\AstroTools-master';
addpath([root '\ASEN5050\tools']);
addpath([root '\ASEN5010\tools']);
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

% And a slightly smaller one to compute a family tangent to feed into the
% SSDR.
fam_tan_alt_diff = 50;
alt_init2 = alt_init - fam_tan_alt_diff; %km
r_init2 = Moon.R+alt_init2;
v_init2 = sqrt(Moon.mu/r_init2);
T_init2 = 2*pi*sqrt(r_init2^3/Moon.mu);

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
        odeset('Events', @y_crossing, 'RelTol', 3e-14, 'AbsTol', 1e-20),propatagor_opts);

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
%     correction = ([STM(1,1) STM(1,5); STM(4,1) STM(4,5)] ...
%         - 1/y_dot*[state_dot(1);state_dot(4)]*[STM(2,1) STM(2,5)])\d;

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

first_PO_state = X_ini;
first_PO_T = T_out(end);

% A second orbit to start a family tangent calculation.

d = [1;1];
tol = 1e-13;
first_guess_state2 = [moon_dist-r_init2;0;0;0;v_init2;0];
X_ini = first_guess_state2;
while abs(d) > tol 
    X = [X_ini; reshape(eye(6),36,1)];

    [T_out,X_out] = ...
        ode45(@SSDR_deriv, [0,P/2], X, ...
        odeset('Events', @y_crossing, 'RelTol', 3e-14, 'AbsTol', 1e-20),propatagor_opts);

    d = -X_out(end,4);
    % STM
    STM = reshape(X_out(end,7:end),6,6);
    y_dot = X_out(end,5);
    state_dot = Lagrange_CR3BP(0,X_out(end,1:6)',propatagor_opts);

    % The correction while holding x constant
    correction = d/( STM(4,5) ...
        - 1/y_dot*state_dot(4)* STM(2,5));
    
    X_ini(5) = X_ini(5) + correction;
    d
end

[Times,X] = ...
    ode45(@Lagrange_CR3BP,[0 10*T_out(end)], X_ini, ...
    ode_opts, propatagor_opts);
plot(X(:,1),X(:,2),'r')

zeroth_PO_state = X_ini;
zeroth_PO_T = T_out(end);

fam_tan_state = first_PO_state-zeroth_PO_state;
fam_tan_T = first_PO_T-zeroth_PO_T;

fam_tan_norm = norm([norm(fam_tan_state) norm(fam_tan_T)]);
fam_tan_norm = norm([fam_tan_state; fam_tan_T]);
fam_tan_state = fam_tan_state/fam_tan_norm;
fam_tan_T = fam_tan_T/fam_tan_norm;

%%
state_in = first_PO_state;
time_in = first_PO_T;
family_tangent = [fam_tan_state;fam_tan_T];
%%
for ii = 1:1000
    ii
    [X_out, T_out] = SSDR(state_in, time_in, family_tangent, 1000, propatagor_opts);
    if mod(ii,50) == 0
        [Times,X] = ...
            ode45(@Lagrange_CR3BP,[0 T_out], X_out, ...
            ode_opts, propatagor_opts);
        % figure
        plot(X(:,1),X(:,2),'k')
    end
    % hold on
    % axis equal
    family_tangent = [X_out-state_in;T_out-time_in]/norm([X_out-state_in;T_out-time_in]);
    state_in = X_out;
    time_in = T_out;
end