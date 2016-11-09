%% Project Integrators
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

% Earth
M1 = Earth.m;

% Moon
M2 = Moon.m;

R = Moon.a;
n = sqrt(G*(M1+M2)/R^3);

propatagor_opts.M1 = M1;
propatagor_opts.M2 = M2;
propatagor_opts.R = R;
propatagor_opts.G = G;
propatagor_opts.nu = M2/(M1+M2);
nu = M2/(M1+M2);
propatagor_opts.n = n;

w = [0 0 n]';
r1 = [-R*M2/(M1+M2) 0 0]';
r2 = [R-R*M2/(M1+M2) 0 0]';
r3 = r1 + [7000 0 0]';
r3 = r2 - [R*(M2/3/M1)^(1/3) 0 0]';

xsi = sym('xsi');
L1_eqn = xsi^5 + (3-nu)*xsi^4 +(3-2*nu)*xsi^3 - nu*xsi^2 -2*nu*xsi - nu;
L1 = solve(L1_eqn);
L1 = double(L1(1));
r3 = r2 - [R*L1 0 0]';


v1 = cross(w,r1);
v2 = cross(w,r2);
v3 = cross(w,r3);

P = 2*pi/n;

ode_opts = odeset('RelTol', 3e-14, 'AbsTol', 1e-20);

%% Test: Numerically integrate primary and secondary bodies
X0 = [r1;v1;r2;v2;r3;v3];
[T,X] = ode45(@Newton3BP,[0 P/4], X0, ode_opts, propatagor_opts);

figure
plot(X(:,1),X(:,2));
hold on
axis equal
plot(X(:,7),X(:,8));
plot(X(:,13),X(:,14));

%% Test: Numerically integrate primary and Third bodies, M2 = 0;
propatagor_opts.M2 = 0;
a_init = 7000;
e_init = 0.01;
i_init = pi/6; %30 degrees
RAAN_init = 0;
w_init = pi/4;
f_init = 0;
[r0, v0 ] = OE2cart( a_init,e_init,i_init,RAAN_init,w_init,f_init,...
    G*M1);
X0 = [r1;[0 0 0]';r2;v2;r0+r1;v0];
[T,X] = ode45(@Newton3BP,[0 P/4], X0, ode_opts, propatagor_opts);

figure
plot3(X(:,1),X(:,2),X(:,3));
hold on
% axis equal
% plot(X(:,7),X(:,8));
plot3(X(:,13),X(:,14),X(:,15),'r');

[rows, cols] = size(X);
OE = zeros(rows,6);
for jj = 1:rows
    [a,e,i,W,w,f] = cart2OE(X(jj,13:15)-X(jj,1:3), X(jj,16:18), G*M1);
    % make w, W continuous
    if w > pi
        w = w - 2*pi;
    end
    if W > pi
        W = W - 2*pi;
    end
    OE(jj,:) = [a,e,i,W,w,f];
end
figure('Position', hw_pub.figPosn);
plot(T,OE(:,1),'LineWidth',hw_pub.lineWidth)
ylim([a-1,a+1]);
ylabel('km','fontsize',hw_pub.fontSize)
xlabel('Time (s)','fontsize',hw_pub.fontSize)
title('a, M_2=0','fontsize',hw_pub.fontSize)
saveas(gcf, ['Figures\' 'OE_Conserved_a'],'epsc')

figure('Position', hw_pub.figPosn);
plot(T,OE(:,2),'LineWidth',hw_pub.lineWidth)
ylim([e-0.001,e+0.001]);
ylabel('--','fontsize',hw_pub.fontSize)
xlabel('Time (s)','fontsize',hw_pub.fontSize)
title('e, M_2=0','fontsize',hw_pub.fontSize)
saveas(gcf, ['Figures\' 'OE_Conserved_e'],'epsc')

figure('Position', hw_pub.figPosn);
plot(T,OE(:,3)*180/pi,'LineWidth',hw_pub.lineWidth)
ylim([i-0.001,i+0.001]*180/pi);
ylabel('Degrees','fontsize',hw_pub.fontSize)
xlabel('Time (s)','fontsize',hw_pub.fontSize)
title('i, M_2=0','fontsize',hw_pub.fontSize)
saveas(gcf, ['Figures\' 'OE_Conserved_i'],'epsc')

figure('Position', hw_pub.figPosn);
plot(T,OE(:,4)*180/pi,'LineWidth',hw_pub.lineWidth)
ylim([W-0.001,W+0.001]*180/pi);
ylabel('Degrees','fontsize',hw_pub.fontSize)
xlabel('Time (s)','fontsize',hw_pub.fontSize)
title('\Omega, M_2=0','fontsize',hw_pub.fontSize)
saveas(gcf, ['Figures\' 'OE_Conserved_W'],'epsc')

figure('Position', hw_pub.figPosn);
plot(T,OE(:,5)*180/pi,'LineWidth',hw_pub.lineWidth)
ylim([w-0.001,w+0.001]*180/pi);
ylabel('Degrees','fontsize',hw_pub.fontSize)
xlabel('Time (s)','fontsize',hw_pub.fontSize)
title('\omega, M_2=0','fontsize',hw_pub.fontSize)
saveas(gcf, ['Figures\' 'OE_Conserved_w'],'epsc')

figure('Position', hw_pub.figPosn);
plot(T,OE(:,6)*180/pi,'LineWidth',hw_pub.lineWidth)
ylabel('Degrees','fontsize',hw_pub.fontSize)
xlabel('Time (s)','fontsize',hw_pub.fontSize)
title('f, M_2=0','fontsize',hw_pub.fontSize)

% percent errors
figure('Position', hw_pub.figPosn);
plot(T,abs((OE(:,1)-OE(1,1))/-OE(1,1))*100,'LineWidth',hw_pub.lineWidth)
ylabel('Percent','fontsize',hw_pub.fontSize)
xlabel('Time (s)','fontsize',hw_pub.fontSize)
title('a percent error, M_2=0','fontsize',hw_pub.fontSize)
saveas(gcf, ['Figures\' 'OE_Conserved_a_PE'],'epsc')

figure('Position', hw_pub.figPosn);
plot(T,abs((OE(:,2)-OE(1,2))/-OE(1,2))*100,'LineWidth',hw_pub.lineWidth)
ylabel('Percent','fontsize',hw_pub.fontSize)
xlabel('Time (s)','fontsize',hw_pub.fontSize)
title('e percent error, M_2=0','fontsize',hw_pub.fontSize)
saveas(gcf, ['Figures\' 'OE_Conserved_e_PE'],'epsc')

figure('Position', hw_pub.figPosn);
plot(T,abs((OE(:,3)-OE(1,3))/-OE(1,3))*100,'LineWidth',hw_pub.lineWidth)
ylabel('Percent','fontsize',hw_pub.fontSize)
xlabel('Time (s)','fontsize',hw_pub.fontSize)
title('i percent error, M_2=0','fontsize',hw_pub.fontSize)
saveas(gcf, ['Figures\' 'OE_Conserved_i_PE'],'epsc')

figure('Position', hw_pub.figPosn);
plot(T,abs((OE(:,4)-OE(1,4)))*100,'LineWidth',hw_pub.lineWidth)
ylabel('Degrees','fontsize',hw_pub.fontSize)
xlabel('Time (s)','fontsize',hw_pub.fontSize)
title('\Omega absolute error, M_2=0','fontsize',hw_pub.fontSize)
saveas(gcf, ['Figures\' 'OE_Conserved_W_PE'],'epsc')

figure('Position', hw_pub.figPosn);
plot(T,abs((OE(:,5)-OE(1,5))/OE(1,5))*100,'LineWidth',hw_pub.lineWidth)
ylabel('Percent','fontsize',hw_pub.fontSize)
xlabel('Time (s)','fontsize',hw_pub.fontSize)
title('\omega percent error, M_2=0','fontsize',hw_pub.fontSize)
saveas(gcf, ['Figures\' 'OE_Conserved_w_PE'],'epsc')


%% Test: Lagrange Vs Newton near L1
propatagor_opts.M2 = M2;
w = [0 0 n]';
r3 = r2 - [R*L1 0 0]';
% r3 = r2 - [R/2 -R*sin(pi/3) 0]'; %L4, which is more stable
v3 = [0 0 0]';
X0 = [r3;v3];

sim_length = P/4;
increment = 3600; %s
int_time = 0:3600:P/4;

[T_Lagrange,X_Lagrange] = ...
    ode45(@Lagrange_CR3BP,int_time, X0, ode_opts, propatagor_opts);
angs = T_Lagrange*n;
X_Lagrange_Inertial = zeros(length(T_Lagrange),3);
Vel_Lagrange_Inertial = zeros(length(T_Lagrange),3);
J = zeros(length(T_Lagrange),1);
for ii = 1:length(T_Lagrange)
    DCM = Euler2DCM('321',[-angs(ii) 0 0]);
    X_Lagrange_Inertial(ii,:) = (DCM*X_Lagrange(ii,1:3)')';
    Vel_Lagrange_Inertial(ii,:) = (DCM*(X_Lagrange(ii,4:6) ...
        + cross(w,X_Lagrange(ii,1:3)))')';
    r_rel = X_Lagrange(ii,1:3);
    v_rel = X_Lagrange(ii,4:6);
    J(ii) = dot(v_rel,v_rel)/2 -dot(cross(w,r_rel),cross(w,r_rel))/2 ...
        -G*M1/norm(r_rel+[nu*R 0 0]) - G*M2/norm(r_rel-[(1-nu)*R 0 0]);
end

v3 = cross(w,r3);
X0 = [r1;v1;r2;v2;r3;v3];
[T,X_Newton] = ode45(@Newton3BP,int_time, X0, ode_opts, propatagor_opts);

figure;
plot(X_Lagrange(:,1), X_Lagrange(:,2))
axis equal

figure;
hold on
axis equal
plot(X_Lagrange_Inertial(:,1), X_Lagrange_Inertial(:,2))
plot(X_Newton(:,13), X_Newton(:,14),'r')
plot(X_Newton(:,7), X_Newton(:,8),'k')

% Error in position
figure('Position', hw_pub.figPosn);
hold on
plot(T_Lagrange,...
    abs((X_Lagrange_Inertial(:,1) - X_Newton(:,13))./X_Newton(:,13))*100,...
    'LineWidth',hw_pub.lineWidth);
plot(T_Lagrange,...
    abs((X_Lagrange_Inertial(:,2) - X_Newton(:,14))./X_Newton(:,14))*100,...
    'r','LineWidth',hw_pub.lineWidth);
plot(T_Lagrange,...
    abs((X_Lagrange_Inertial(:,3) - X_Newton(:,15))./X_Newton(:,15))*100,...
    'k','LineWidth',hw_pub.lineWidth);
title('Error in Inertial Position','fontsize',hw_pub.fontSize)
ylabel('Percent Error','fontsize',hw_pub.fontSize)
xlabel('Time','fontsize',hw_pub.fontSize)
legend('x','y','z','fontsize',hw_pub.fontSize)
saveas(gcf, ['Figures\' 'Pos_Err'],'epsc')

% Error in Velocity
figure('Position', hw_pub.figPosn);
hold on
plot(T_Lagrange,...
    abs((Vel_Lagrange_Inertial(:,1) - X_Newton(:,16))./X_Newton(:,16))*100,...
    'LineWidth',hw_pub.lineWidth);
plot(T_Lagrange,...
    abs((Vel_Lagrange_Inertial(:,2) - X_Newton(:,17))./X_Newton(:,17))*100,...
    'r','LineWidth',hw_pub.lineWidth);
plot(T_Lagrange,...
    abs((Vel_Lagrange_Inertial(:,3) - X_Newton(:,18))./X_Newton(:,18))*100,...
    'k','LineWidth',hw_pub.lineWidth);
title('Error in Inertial Velocity','fontsize',hw_pub.fontSize)
ylabel('Percent Error','fontsize',hw_pub.fontSize)
xlabel('Time','fontsize',hw_pub.fontSize)
legend({'$\dot{x}$','$\dot{y}$','$\dot{z}$'},'Interpreter','latex','fontsize',hw_pub.fontSize)
saveas(gcf, ['Figures\' 'Vel_Err'],'epsc')

% Check that J is conserved
% figure
% plot(J)
figure('Position', hw_pub.figPosn);
plot(T_Lagrange,abs((J-J(1))/J(1))*100,'LineWidth',hw_pub.lineWidth);
ylabel('J Percent Error','fontsize',hw_pub.fontSize);
xlabel('Time, s','fontsize',hw_pub.fontSize);
saveas(gcf, ['Figures\' 'J_conserved'],'epsc')

%% Test: Use the previous earth orbit and verify errors, OE changes.

w = [0 0 n]';
propatagor_opts.M2 = M2;
sim_length = P/4;
increment = 3600; %s
int_time = 0:3600:sim_length;

X0 = [r1;v1;r2;v2;r0+r1;v0+v1];
[T,X_Newton] = ode45(@Newton3BP,int_time, X0, ode_opts, propatagor_opts);

X0 = [r0+r1;v0+v1-cross(w,r0+r1)];
[T_Lagrange,X_Lagrange] = ...
    ode45(@Lagrange_CR3BP,int_time, X0, ode_opts, propatagor_opts);

angs = T_Lagrange*n;
X_Lagrange_Inertial = zeros(length(T_Lagrange),3);
Vel_Lagrange_Inertial = zeros(length(T_Lagrange),3);
J = zeros(length(T_Lagrange),1);
for ii = 1:length(T_Lagrange)
    DCM = Euler2DCM('321',[-angs(ii) 0 0]);
    X_Lagrange_Inertial(ii,:) = (DCM*X_Lagrange(ii,1:3)')';
    Vel_Lagrange_Inertial(ii,:) = (DCM*(X_Lagrange(ii,4:6) ...
        + cross(w,X_Lagrange(ii,1:3)))')';
    r_rel = X_Lagrange(ii,1:3);
    v_rel = X_Lagrange(ii,4:6);
    J(ii) = dot(v_rel,v_rel)/2 -dot(cross(w,r_rel),cross(w,r_rel))/2 ...
        -G*M1/norm(r_rel+[nu*R 0 0]) - G*M2/norm(r_rel-[(1-nu)*R 0 0]);
end

figure;
hold on
axis equal
plot(X_Lagrange_Inertial(:,1), X_Lagrange_Inertial(:,2))
plot(X_Newton(:,13), X_Newton(:,14),'r')
% plot(X_Newton(:,7), X_Newton(:,8),'k')

% Error in position
figure('Position', hw_pub.figPosn);
hold on
plot(T_Lagrange,...
    abs((X_Lagrange_Inertial(:,1) - X_Newton(:,13))./X_Newton(:,13))*100,...
    'LineWidth',hw_pub.lineWidth);
plot(T_Lagrange,...
    abs((X_Lagrange_Inertial(:,2) - X_Newton(:,14))./X_Newton(:,14))*100,...
    'r','LineWidth',hw_pub.lineWidth);
plot(T_Lagrange,...
    abs((X_Lagrange_Inertial(:,3) - X_Newton(:,15))./X_Newton(:,15))*100,...
    'k','LineWidth',hw_pub.lineWidth);
title('Error in Inertial Position','fontsize',hw_pub.fontSize)
ylabel('Percent Error','fontsize',hw_pub.fontSize)
xlabel('Time','fontsize',hw_pub.fontSize)
legend('x','y','z','fontsize',hw_pub.fontSize)
saveas(gcf, ['Figures\' 'Pos_Err_Example_Orb'],'epsc')

% Error in Velocity
figure('Position', hw_pub.figPosn);
hold on
plot(T_Lagrange,...
    abs((Vel_Lagrange_Inertial(:,1) - X_Newton(:,16))./X_Newton(:,16))*100,...
    'LineWidth',hw_pub.lineWidth);
plot(T_Lagrange,...
    abs((Vel_Lagrange_Inertial(:,2) - X_Newton(:,17))./X_Newton(:,17))*100,...
    'r','LineWidth',hw_pub.lineWidth);
plot(T_Lagrange,...
    abs((Vel_Lagrange_Inertial(:,3) - X_Newton(:,18))./X_Newton(:,18))*100,...
    'k','LineWidth',hw_pub.lineWidth);
title('Error in Inertial Velocity','fontsize',hw_pub.fontSize)
ylabel('Percent Error','fontsize',hw_pub.fontSize)
xlabel('Time','fontsize',hw_pub.fontSize)
legend({'$\dot{x}$','$\dot{y}$','$\dot{z}$'},'Interpreter','latex','fontsize',hw_pub.fontSize)
saveas(gcf, ['Figures\' 'Vel_Err_Example_Orb'],'epsc')

% Check that J is conserved
figure('Position', hw_pub.figPosn);
plot(T_Lagrange,abs((J-J(1))/J(1))*100,'LineWidth',hw_pub.lineWidth);
ylabel('J Percent Error','fontsize',hw_pub.fontSize);
xlabel('Time, s','fontsize',hw_pub.fontSize);

[rows, cols] = size(X_Newton);
OE = zeros(rows,6);
for jj = 1:rows
    [a,e,i,W,w,f] = cart2OE(X_Newton(jj,13:15)-X_Newton(jj,1:3), X_Newton(jj,16:18), G*M1);
    % make w, W continuous
    if w > pi
        w = w - 2*pi;
    end
    if W > pi
        W = W - 2*pi;
    end
    OE(jj,:) = [a,e,i,W,w,f];
end
figure('Position', hw_pub.figPosn);
plot(T,OE(:,1),'LineWidth',hw_pub.lineWidth)
ylabel('km','fontsize',hw_pub.fontSize)
xlabel('Time (s)','fontsize',hw_pub.fontSize)
title('a, M_2=0','fontsize',hw_pub.fontSize)