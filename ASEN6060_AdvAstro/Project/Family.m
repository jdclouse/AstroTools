%% Project Find Family
% John Clouse
close all
clear all

% Initialize some things.
root = 'C:\Users\John\Documents\Astro';
% root = 'C:\Users\jclouse1\Documents\Class\AstroTools-master';
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

nu = M2/(M1+M2);
propatagor_opts.nu = M2/(M1+M2);
propatagor_opts.M1 = M1;
propatagor_opts.M2 = M2;
propatagor_opts.R = R;
propatagor_opts.G = G;
propatagor_opts.n = n;
propatagor_opts.M1 = 1-nu;
propatagor_opts.M2 = nu;
propatagor_opts.R = 1;
propatagor_opts.G = 1;
propatagor_opts.n = 1;

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
% Canonical units
first_guess_state = first_guess_state/R;
first_guess_state(4:6) = first_guess_state(4:6)*P/(2*pi);
int_time = [0 1/2];
[Times,X] = ...
    ode45(@Lagrange_CR3BP,int_time, first_guess_state, ...
    ode_opts, propatagor_opts);

figure
plot(X(:,1),X(:,2))
hold on
axis equal
% figure
% plot(X(:,3))
% hold on
% axis equal
d = [1;1];
% tol = 1e-13; %worked for earth-moon non-canonical
tol = 1e-12;
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
% tol = 1e-13;
first_guess_state2 = [moon_dist-r_init2;0;0;0;v_init2;0];
% Canonical units
first_guess_state2 = first_guess_state2/R;
first_guess_state2(4:6) = first_guess_state2(4:6)*P/(2*pi);
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
for ii = 1:480
    ii
    [X_out, T_out] = SSDR(state_in, time_in, family_tangent, .01, propatagor_opts);
    if mod(ii,20) == 0
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

%% Stability
X_ini2 = state_in;
X_ini2(2) = 0;
X_ini2(3) = 0;
X_ini2(4) = 0;
X_ini2(6) = 0;
d = 1;
tol = 1e-14
while abs(d) > tol 
    X = [X_ini2; reshape(eye(6),36,1)];

    [T_out,X_out] = ...
        ode45(@SSDR_deriv, [0,P/2], X, ...
        odeset('Events', @y_crossing, ...
        'RelTol', 3e-14, 'AbsTol', 1e-20),propatagor_opts);

    d = -X_out(end,4);
    % STM
    STM = reshape(X_out(end,7:end),6,6);
    y_dot = X_out(end,5);
    state_dot = Lagrange_CR3BP(0,X_out(end,1:6)',propatagor_opts);

    % The correction while holding x constant
    correction = d/( STM(4,5) ...
        - 1/y_dot*state_dot(4)* STM(2,5));
    
    X_ini2(5) = X_ini2(5) + correction;
    d
end
stab_T = T_out(end);
stab_state_i = X_ini2;
stab_state_i(5) = stab_state_i(5)+ 0.01*P/(2*pi)/R;
% [T_stab,X_stab] = ode45(@SSDR_deriv,[0 20*stab_T], ...
[T_stab,X_stab] = ode45(@SSDR_deriv,[0 12*20*2*pi], ...
    [stab_state_i; reshape(eye(6),36,1)], ...
    ode_opts, propatagor_opts);
stab_STM = reshape(X_stab(end,7:end),6,6);
e_vals = eigs(stab_STM);

figure
plot(X_stab(:,1),X_stab(:,2))

DRO_stab_plot = figure;
ang_circ = linspace(0,2*pi,100);
x_circ = cos(ang_circ);
y_circ = sin(ang_circ);
plot(x_circ, y_circ)
axis equal
hold on
for ii = 1:length(e_vals)
    plot(real(e_vals(ii)), imag(e_vals(ii)),'rx')
end
%%
dyn_mat = hamiltonian_dyn_mat(stab_state_i, propatagor_opts);
dyn_mat_e_vals = eigs(dyn_mat)

%% Prograde
d = 1;
tol = 1e-13;
X_ini = first_guess_state;
X_ini(5) = -X_ini(5);
while abs(d) > tol 
    X = [X_ini; reshape(eye(6),36,1)];

    [T_out,X_out] = ...
        ode45(@SSDR_deriv, [0,P/2], X, ...
        odeset('Events', @neg_y_crossing, 'RelTol', 3e-14, 'AbsTol', 1e-20),propatagor_opts);

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
figure
plot(X(:,1),X(:,2))
hold on
axis equal

first_DPO_state = X_ini;
first_DPO_T = T_out(end);

% A second orbit to start a family tangent calculation.

d = 1
X_ini = first_guess_state2;
X_ini(5) = -X_ini(5);
while abs(d) > tol 
    X = [X_ini; reshape(eye(6),36,1)];

    [T_out,X_out] = ...
        ode45(@SSDR_deriv, [0,P/2], X, ...
        odeset('Events', @neg_y_crossing, 'RelTol', 3e-14, 'AbsTol', 1e-20),propatagor_opts);

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

zeroth_DPO_state = X_ini;
zeroth_DPO_T = T_out(end);

DPO_fam_tan_state = first_DPO_state-zeroth_DPO_state;
DPO_fam_tan_T = first_DPO_T-zeroth_DPO_T;

DPO_fam_tan_norm = norm([DPO_fam_tan_state; DPO_fam_tan_T]);
DPO_fam_tan_state = DPO_fam_tan_state/DPO_fam_tan_norm;
DPO_fam_tan_T = DPO_fam_tan_T/DPO_fam_tan_norm;
%% Initialize prograde states
DPO_state_in = first_DPO_state;
DPO_time_in = first_DPO_T;
DPO_family_tangent = [DPO_fam_tan_state;DPO_fam_tan_T];
%%
for ii = 1:300
    fprintf('%d\n',ii)
    [X_out, T_out] = SSDR(DPO_state_in, DPO_time_in, DPO_family_tangent, .01, propatagor_opts);
    if mod(ii,20) == 0
        [Times,X] = ...
            ode45(@Lagrange_CR3BP,[0 T_out], X_out, ...
            ode_opts, propatagor_opts);
        % figure
        plot(X(:,1),X(:,2),'k')
        drawnow;
    end
    % hold on
    % axis equal
    DPO_family_tangent = [X_out-DPO_state_in;T_out-DPO_time_in]/norm([X_out-DPO_state_in;T_out-DPO_time_in]);
    DPO_state_in = X_out;
    DPO_time_in = T_out;
end

%%
prograde_state = DPO_state_in;
prograde_state(5) = prograde_state(5) + 0.01*P/(2*pi)/R;
% [T_stab,X_pro_stab] = ode45(@SSDR_deriv,[0 20*DPO_time_in], ...
[T_stab,X_pro_stab] = ode45(@SSDR_deriv,[0 12*20*2*pi], ...
    [prograde_state; reshape(eye(6),36,1)], ...
    ode_opts, propatagor_opts);
prog_stab_STM = reshape(X_pro_stab(end,7:end),6,6);
e_vals_prog = eigs(prog_stab_STM);
DPO_stab_plot = figure;
plot(x_circ, y_circ)
axis equal
hold on
for ii = 1:length(e_vals_prog)
    plot(real(e_vals_prog(ii)), imag(e_vals_prog(ii)),'rx')
end
prograde_orb_plot = figure;
plot(X_pro_stab(:,1),X_pro_stab(:,2))
axis equal


