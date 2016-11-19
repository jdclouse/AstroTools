%% ASEN 5227 Project
% John Clouse
clear all
close all

scenarios = {'Part 1', 'Part 2'};

for scenario = scenarios

if strcmp(scenario, 'Part 1')
    num_pts = 1000;
    time = linspace(0,2*pi,num_pts);
else
    mav_data = dlmread('mav_data_txt.txt');
    time = mav_data(:,4)';
    num_pts = length(time);
end
t_diffs = (time(2:end) - time(1:end-1));
R = [time.*cos(time); time.*sin(2*time); time];
V = [(cos(time)-time.*sin(time)); ...
    (sin(2*time)+2*time.*cos(2*time));...
    1*ones(1,num_pts)];
A = [(-2*sin(time)-time.*cos(time));...
    (4*cos(2*time)-4*time.*sin(2*time));...
    zeros(1,num_pts)];

if strcmp(scenario, 'Part 1')
    r_bar = [cos(time); sin(2*time); cos(2*time)]; %Frenet
    v_bar = [-sin(time); 2*cos(2*time); -2*sin(2*time)];
    a_bar = [-cos(time); -4*sin(2*time); -4*cos(2*time)];
else
    r_bar = mav_data(:,1:3)';
end

% Forward differentiation
% Mothership velocity
V_forward = forward_diff(R,time);
figure;
hold on
plot(abs(V(1,1:end-1) - V_forward(1,:)));
plot(abs(V(2,1:end-1) - V_forward(2,:)),'r');
plot(abs(V(3,1:end-1) - V_forward(3,:)),'k');
title('MS Velocity: Error Between Analytical and Numerical Solution')

% Mothership acceleration
A_forward = forward_diff(V_forward,time);
figure;
hold on
plot(abs(A(1,1:end-2) - A_forward(1,:)));
plot(abs(A(2,1:end-2) - A_forward(2,:)),'r');
plot(abs(A(3,1:end-2) - A_forward(3,:)),'k');
title('MS Acceleration: Error Between Analytical and Numerical Solution')

v_bar_forward = forward_diff(r_bar,time);
a_bar_forward = forward_diff(v_bar_forward,time);
if strcmp(scenario, 'Part 1')
    % MAV velocity
    figure;
    hold on
    plot(abs(v_bar(1,1:end-1) - v_bar_forward(1,:)));
    plot(abs(v_bar(2,1:end-1) - v_bar_forward(2,:)),'r');
    plot(abs(v_bar(3,1:end-1) - v_bar_forward(3,:)),'k');
    title('MAV Velocity: Error Between Analytical and Numerical Solution')

    % MAV acceleration
    figure;
    hold on
    plot(abs(a_bar(1,1:end-2) - a_bar_forward(1,:)));
    plot(abs(a_bar(2,1:end-2) - a_bar_forward(2,:)),'r');
    plot(abs(a_bar(3,1:end-2) - a_bar_forward(3,:)),'k');
    title('MAV Acceleration: Error Between Analytical and Numerical Solution')
end

% Mothership's local Frenet frame
t = zeros(3,num_pts);
b = zeros(3,num_pts);
n = zeros(3,num_pts);
for ii = 1:num_pts
    t(:,ii) = V(:,ii)/norm(V(:,ii));
    b(:,ii) = cross(V(:,ii),A(:,ii))/norm(cross(V(:,ii),A(:,ii)));
    n(:,ii) = cross(b(:,ii),t(:,ii));
end

if strcmp(scenario, 'Part 1')
    figure;
    plot3(R(1,:),R(2,:),R(3,:));
    xlabel('i'); ylabel('j'); zlabel('k')
    hold on
    plot_idx = [1:20:num_pts num_pts];
    quiver3(R(1,plot_idx),R(2,plot_idx),R(3,plot_idx),...
        t(1,plot_idx),t(2,plot_idx),t(3,plot_idx),'r');
    quiver3(R(1,plot_idx),R(2,plot_idx),R(3,plot_idx),...
        b(1,plot_idx),b(2,plot_idx),b(3,plot_idx),'r');
    quiver3(R(1,plot_idx),R(2,plot_idx),R(3,plot_idx),...
        n(1,plot_idx),n(2,plot_idx),n(3,plot_idx),'r');
    xlabel('i'); ylabel('j'); zlabel('k')
    title('MS observed from GS with Frenet frame vectors')
end

% MAV's local Frenet frame
t_MAV = zeros(3,num_pts);
b_MAV = zeros(3,num_pts);
n_MAV = zeros(3,num_pts);
for ii = 1:num_pts-2 % missing 2 points at end from 2 forward diffs
    t_MAV(:,ii) = v_bar_forward(:,ii)/norm(v_bar_forward(:,ii));
    b_MAV(:,ii) = cross(v_bar_forward(:,ii),a_bar_forward(:,ii))...
        /norm(cross(v_bar_forward(:,ii),a_bar_forward(:,ii)));
    n_MAV(:,ii) = cross(b_MAV(:,ii),t_MAV(:,ii));
end


if strcmp(scenario, 'Part 1')
    plot_idx = [1:20:num_pts num_pts];
else
    plot_idx=[1:4:num_pts num_pts];
end
figure;
plot3(r_bar(1,:),r_bar(2,:),r_bar(3,:));
xlabel('i'); ylabel('j'); zlabel('k')
hold on
quiver3(r_bar(1,plot_idx),r_bar(2,plot_idx),r_bar(3,plot_idx),...
    t_MAV(1,plot_idx),t_MAV(2,plot_idx),t_MAV(3,plot_idx),'r');
quiver3(r_bar(1,plot_idx),r_bar(2,plot_idx),r_bar(3,plot_idx),...
    b_MAV(1,plot_idx),b_MAV(2,plot_idx),b_MAV(3,plot_idx),'r');
quiver3(r_bar(1,plot_idx),r_bar(2,plot_idx),r_bar(3,plot_idx),...
    n_MAV(1,plot_idx),n_MAV(2,plot_idx),n_MAV(3,plot_idx),'r');
xlabel('t'); ylabel('n'); zlabel('b')
title('MAV observed from GS with Frenet frame vectors')


if strcmp(scenario, 'Part 1')
    % Mothership normal, tangent accels
    MS_accel_tangent = zeros(1,length(A_forward));
    MS_accel_normal = zeros(1,length(A_forward));
    MS_accel_bi = zeros(1,length(A_forward));
    for ii = 1:length(A_forward)
        MS_accel_tangent(ii) = dot(t(:,ii),A_forward(:,ii));
        MS_accel_normal(ii) = dot(n(:,ii),A_forward(:,ii));
        MS_accel_bi(ii) = dot(b(:,ii),A_forward(:,ii));
    end

    figure;
    plot(MS_accel_tangent)
    title('MS Tangential Acceleration wrt GS')
    figure;
    plot(MS_accel_normal)
    title('MS Normal Acceleration wrt GS')
    figure;
    plot(MS_accel_bi)
title('MS Binormal Acceleration wrt GS')
end

% MAV normal, tangent accels
MAV_accel_tangent = zeros(1,length(a_bar_forward));
MAV_accel_normal = zeros(1,length(a_bar_forward));
MAV_accel_bi = zeros(1,length(a_bar_forward));
for ii = 1:length(a_bar_forward)
    MAV_accel_tangent(ii) = dot(t_MAV(:,ii),a_bar_forward(:,ii));
    MAV_accel_normal(ii) = dot(n_MAV(:,ii),a_bar_forward(:,ii));
    MAV_accel_bi(ii) = dot(b_MAV(:,ii),a_bar_forward(:,ii));
end

figure;
plot(MAV_accel_tangent)
title('MAV Tangential Acceleration wrt MS')
figure;
plot(MAV_accel_normal)
title('MAV Normal Acceleration wrt MS')
figure;
plot(MAV_accel_bi)
title('MAV Binormal Acceleration wrt MS')

% MAV wrt GS
% The 
r = nan(3,length(t)); % nans so that the zeros don't drag down the plots
v = zeros(3,length(t));
phi = zeros(1,length(t));
theta = zeros(1,length(t));
psi = zeros(1,length(t));
phi_dot = zeros(1,length(t));
theta_dot = zeros(1,length(t));
psi_dot = zeros(1,length(t));
w_MS_Frenet_wrt_GS = zeros(3,length(t));
w_MS_Frenet_wrt_body = zeros(3,length(t));
for ii = 1:length(t) - 2 % due to forward diff
    rot = [t(:,ii)';n(:,ii)';b(:,ii)'];
    % Orthogonal transformation, so the inverse is the transpose.
    G = rot';
    r(:,ii) = R(:,ii)+G*r_bar(:,ii);
    phi(ii) = atan2(G(1,3),-G(2,3));
    theta(ii) = atan2(sqrt(1-G(3,3)^2),G(3,3));
    psi(ii) = atan2(G(3,1),G(3,2));
    if phi(ii) < 0
        phi(ii) = phi(ii) + 2*pi;
    end
    if theta(ii) < 0
        theta(ii) = theta(ii) + 2*pi;
    end
    if psi(ii) < 0
        psi(ii) = psi(ii) + 2*pi;
    end
    if ii >= 2 && abs(phi(ii)-phi(ii-1)) > pi
        phi(ii) = phi(ii) + 2*pi;
    end
    if ii >= 2 && abs(theta(ii)-theta(ii-1)) > pi
        theta(ii) = theta(ii) + 2*pi;
    end
    if ii >= 2 && abs(psi(ii)-psi(ii-1)) > pi
        psi(ii) = psi(ii) + 2*pi;
    end
    if ii >=2 % for first-order backward differencing
        phi_dot(ii) = (phi(ii)-phi(ii-1))/t_diffs(ii-1);
        theta_dot(ii) = (theta(ii)-theta(ii-1))/t_diffs(ii-1);
        psi_dot(ii) = (psi(ii)-psi(ii-1))/t_diffs(ii-1);
    end
    w_MS_Frenet_wrt_GS(1,ii) = psi_dot(ii)*sin(theta(ii))*sin(phi(ii)) ...
        + theta_dot(ii)*cos(phi(ii));
    w_MS_Frenet_wrt_GS(2,ii) = -psi_dot(ii)*sin(theta(ii))*cos(phi(ii)) ...
        + theta_dot(ii)*sin(phi(ii));
    w_MS_Frenet_wrt_GS(3,ii) = psi_dot(ii)*cos(theta(ii))+phi_dot(ii);
    w_MS_Frenet_wrt_body(1,ii) = phi_dot(ii)*sin(theta(ii))*sin(psi(ii))...
        +theta_dot(ii)*cos(psi(ii));
    w_MS_Frenet_wrt_body(2,ii) = phi_dot(ii)*sin(theta(ii))*cos(psi(ii))...
        -theta_dot(ii)*sin(psi(ii));
    w_MS_Frenet_wrt_body(3,ii) = phi_dot(ii)*cos(theta(ii))+psi_dot(ii);
    
    v(:,ii) = V(:,ii) + G*(v_bar_forward(:,ii) ...
        + cross(w_MS_Frenet_wrt_body(:,ii),r_bar(:,ii)));
end
figure
plot3(r(1,:),r(2,:),r(3,:));
hold on
plot3(R(1,:),R(2,:),R(3,:),'r');
xlabel('i'); ylabel('j'); zlabel('k')
title('MAV observed from GS')

v_forward = forward_diff(r,time);
figure
hold on
plot(abs(v(1,1:end-1) - v_forward(1,:)));
plot(abs(v(2,1:end-1) - v_forward(2,:)),'r');
plot(abs(v(3,1:end-1) - v_forward(3,:)),'k');
title('MAV Velocity wrt GS: Error Between Euler-Angle-Rotation and Numerical Differentiation Solution')
end