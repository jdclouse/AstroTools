%% P-Map

%%
% find an orbit at a similar distance for ARM:
% ARM_r = (moon_dist-(70000+Moon.R))/R;
% for ii = 2:store_cnt-1
%     if X_PO_store(1,ii) < ARM_r
%         member_PO_prograde = X_PO_store(:,ii-1);
%         member_T_prograde = T_PO_store(ii-1);
%         break;
%     end
% end

member_PO_prograde = X_PO_store(:,13);
member_T_prograde = T_PO_store(13);

%% Perturb X
x_only_p_map_x_vx_prograde = figure(Position', hw_pub.figPosn);
x_only_p_map_prograde = figure(Position', hw_pub.figPosn);
for ii = 4
    tic
x_only_pert_prograde = member_PO_prograde + disturbances(ii)/R*[1 0 0 0 0 0 ]';
x_only_pert_prograde(5) = find_ydot(member_PO_prograde,x_only_pert_prograde);
[~,~,~,X_e_xo_prograde,~] = ode45(@Lagrange_CR3BP,[0 member_T_prograde*1000], x_only_pert_prograde, ...
    odeset('Events', @y_crossing_1k, 'RelTol', 3e-14, 'AbsTol', 1e-14), ...
    propatagor_opts);
% C = compute_C(x_only_pert(1:3),x_only_pert(4:6))

% plot(X_e_xo(:,2))
figure(x_only_p_map_x_vx_prograde);
plot(X_e_xo_prograde(:,1),X_e_xo_prograde(:,3),'x','color', colors{ii})
axis equal
hold on
drawnow
figure(x_only_p_map_prograde);
plot(X_e_xo_prograde(:,1),X_e_xo_prograde(:,4),'x','color', colors{ii})
axis equal
hold on
drawnow
toc
% plot(X_e(:,1),X_e(:,6),'x')
end

figure
plot(X_e_xo_prograde(:,1),X_e_xo_prograde(:,2))


%% Perturb Z
Z_only_p_map_prograde = figure(Position', hw_pub.figPosn);
z_only_p_map_z_vz_prograde = figure(Position', hw_pub.figPosn);
for ii = 1:length(disturbances)
    tic
Z_only_pert_prograde = member_PO_prograde + disturbances(ii)/R*[0 0 1 0 0 0 ]';
Z_only_pert_prograde(5) = find_ydot(member_PO_prograde,Z_only_pert_prograde);
[~,X_zo_prograde,~,X_e_zo_prograde,~] = ode45(@Lagrange_CR3BP,[0 member_T_prograde*1000], Z_only_pert_prograde, ...
    odeset('Events', @y_crossing_1k, 'RelTol', 3e-14, 'AbsTol', 1e-14), ...
    propatagor_opts);

figure(Z_only_p_map_prograde);
plot(X_e_zo_prograde(:,1),X_e_zo_prograde(:,3),'x','color', colors{ii})
axis equal; hold on; drawnow;
figure(z_only_p_map_z_vz_prograde);
plot(X_e_zo_prograde(:,3),X_e_zo_prograde(:,6),'x','color', colors{ii})
axis equal; hold on; drawnow;
toc
if disturbances(ii) == 1000
    figure(Position', hw_pub.figPosn);
    plot3(X_zo_prograde(:,1),X_zo_prograde(:,2),X_zo_prograde(:,3))
end
end
%% Perturb VX
vx_only_pert_prograde = member_PO_prograde + 0.01/R*P/(2*pi)*[0 0 0 1 0 0 ]';
vx_only_pert_prograde(5) = find_ydot(member_PO_prograde,vx_only_pert_prograde);
[~,~,~,X_e_vxo_prograde,~] = ode45(@Lagrange_CR3BP,[0 member_T_prograde*1000], vx_only_pert_prograde, ...
    odeset('Events', @y_crossing_1k, 'RelTol', 3e-14, 'AbsTol', 1e-14), ...
    propatagor_opts);

vx_only_p_map_prograde = figure(Position', hw_pub.figPosn);
plot(X_e_vxo_prograde(:,1),X_e_vxo_prograde(:,3),'x')
axis equal
vx_only_p_map_x_vx_prograde = figure(Position', hw_pub.figPosn);
plot(X_e_vxo_prograde(:,1),X_e_vxo_prograde(:,4),'x')
axis equal
% plot(X_e(:,1),X_e(:,6),'x')
%% Perturb XZ
xz_pert_prograde = member_PO_prograde + 1000/R/sqrt(2)*[1 0 1 0 0 0 ]';
xz_pert_prograde(5) = find_ydot(member_PO_prograde,xz_pert_prograde);
[~,~,~,X_e_xz_prograde,~] = ode45(@Lagrange_CR3BP,[0 member_T_prograde*1000], xz_pert_prograde, ...
    odeset('Events', @y_crossing_1k, 'RelTol', 3e-14, 'AbsTol', 1e-14), ...
    propatagor_opts);

xz_p_map_prograde = figure(Position', hw_pub.figPosn);
% plot(X_e_xo(:,2))
plot(X_e_xz_prograde(:,1),X_e_xz_prograde(:,3),'x')
axis equal
xz_p_map_x_vx_prograde = figure(Position', hw_pub.figPosn);
plot(X_e_xz_prograde(:,1),X_e_xz_prograde(:,4),'x')
axis equal