%% B-plane target
ode_opts = odeset('RelTol', 1e-12, 'AbsTol', 1e-20,...
    'Events',@stop_int);

output_state = EKFoutput.state_store(:,end);
final_P = EKFoutput.final_P;

% Get the integration times
processed_obs = ObsMassaged(obs_to_process,:);

[T, X_to_SOI] = ode45(@flyby_two_body_state_dot, ...
    [processed_obs(end,2), processed_obs(end,2)+400*86400], ...
    [output_state; reshape(eye(7),49,1)], ...
        ode_opts, filter_opts.propagator_opts);
    
v_inf = norm(X_to_SOI(end,4:6));

f = acosh(1+ v_inf*v_inf/mu_earth* norm(X_to_SOI(end,1:3)));

LTOF = mu_earth/v_inf/v_inf/v_inf*(sinh(f)-f);

S_hat = X_to_SOI(end,4:6)/v_inf;

% T_hat
% The B Plane. note that without a V_out, this isn't very correct except
% the B Plane definition (so keep only that).
[~, ~, B_plane, ~, ~] = ...
    BPlaneTarget(X_to_SOI(end,4:6)', X_to_SOI(end,4:6)', mu_earth);

ode_opts = odeset('RelTol', 1e-12, 'AbsTol', 1e-20);
[~, X_to_BPlane] = ode45(@flyby_two_body_state_dot, ...
    [processed_obs(end,2), T(end)+LTOF], ...
    [output_state; reshape(eye(7),49,1)], ...
        ode_opts, filter_opts.propagator_opts);

% aim point: r_inf projected on the B plane
BT = dot(X_to_SOI(end,1:3)',B_plane(:,2));
BR = dot(X_to_SOI(end,1:3)',B_plane(:,3));

STM_DCO_intercept = reshape(X_to_BPlane(end,8:end)',7,7);
P_intercept = STM_DCO_intercept*final_P*STM_DCO_intercept';
P_int_Bplane = B_plane*P_intercept(1:3,1:3)*B_plane';
[U,D] = eig(P_int_Bplane(2:3,2:3));
ellip_tilt = atan2(U(2,1),-U(1,1));

a=3*sqrt(max(max(D))); % horizontal radius
b=3*sqrt(min(max(D))); % vertical radius
x0=BT; % x0,y0 ellipse centre coordinates
y0=BR;
t=-pi:0.01:pi;
x=a*cos(t);
y=b*sin(t);
ell_rot = [cos(ellip_tilt) -sin(ellip_tilt); sin(ellip_tilt) cos(ellip_tilt)];
coords_prime = [];
for ii = 1:length(t)
    coords_prime(:,ii) = ell_rot*[x(ii); y(ii)];
end

figure(ell_plot);
hold on
% plot(x0+coords_prime(1,:),y0+coords_prime(2,:),'Color',color_order(kk,:),...
%     'LineWidth',hw_pub.lineWidth)
% plot_handles1 = [plot_handles1 plot(x0, y0, 'x','Color',color_order(kk,:),...
%     'LineWidth',hw_pub.lineWidth)];
plot(x0+coords_prime(1,:),y0+coords_prime(2,:))
plot_handles1 = [plot_handles1 plot(x0, y0, 'x')];
axis equal

figure(ell_plot_co);
hold on

% plot_handles2 = [plot_handles2 plot(coords_prime(1,:),coords_prime(2,:),...
%     'Color',color_order(kk,:),...
%     'LineWidth',hw_pub.lineWidth)];
% axis equal