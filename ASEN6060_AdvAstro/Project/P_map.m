%% P-Map

%%
% find an orbit at a similar distance for ARM:
ARM_r = (moon_dist-(70000+Moon.R))/R;
for ii = 2:store_cnt-1
    if X_PO_store(1,ii) < ARM_r
        member_PO = X_PO_store(:,ii-1);
        member_T = T_PO_store(ii-1);
        member_STM = STM_PO_store(:,:,ii-1);
        break;
    end
end

% a quick integration of the main solution...
[~,X_Lagrange] = ode45(@SSDR_deriv,[0 member_T], [member_PO; reshape(eye(6),36,1)], ...
    ode_opts, propatagor_opts);
C = zeros(1,length(X_Lagrange));
w = [0 0 1]';
compute_U = @(r_rel) propatagor_opts.G*propatagor_opts.M1/norm(r_rel+[nu 0 0]')...
    + propatagor_opts.G*propatagor_opts.M2/norm(r_rel-[(1-nu) 0 0]');
compute_C = @(r_rel,v_rel) dot(v_rel,v_rel)/2  -dot(cross(w,r_rel),cross(w,r_rel))/2 ...
    - compute_U(r_rel);
for ii = 1:length(C)
r_rel = X_Lagrange(ii,1:3)';
v_rel = X_Lagrange(ii,4:6)';
C(ii) = compute_C(r_rel,v_rel);
end
figure
plot(abs(C-C(1))/-C(1))

member_PO_STM = reshape(X_Lagrange(end,7:end),6,6);
member_PO_eig_vals = eigs(member_PO_STM);
[member_PO_eig_vecx,~] = eigs(member_PO_STM);

freq1 = abs(imag(log(member_PO_eig_vals(2))));
freq2 = abs(imag(log(member_PO_eig_vals(4))));

figure; 
plot(X_Lagrange(:,1),X_Lagrange(:,2)); 
axis equal; hold on;
% plot(x_circ*Earth.R/R-nu,y_circ*Earth.R/R, 'k');
fill(x_circ*Earth.R/R-nu,y_circ*Earth.R/R, 'k');
fill(x_circ*Moon.R/R+1-nu,y_circ*Moon.R/R, 'r');
plot(x_circ*Moon.R/R+1-nu,y_circ*Moon.R/R, 'r', 'linewidth', 1.5);
xlim([-0.05 1.2])

%%
find_ydot = @(ref_PO,pert_PO) ...
    sqrt(2*(compute_C(ref_PO(1:3),ref_PO(4:6))...
    +compute_U(pert_PO(1:3)))...
    - pert_PO(4)^2-pert_PO(6)^2 + pert_PO(1)^2);

%% Set up some perturbations, colors
disturbances = [1, 10, 100, 1000, 10000];
colors = {[19 19 23] ./ 255, [235 19 51] ./ 255, [12 130 62] ./ 255, ...
    [72 30 104] ./ 255, [64 99 196] ./ 255};
%% Perturb X
x_only_p_map_x_z = figure(Position', hw_pub.figPosn);
x_only_p_map_x_vx = figure(Position', hw_pub.figPosn);
for ii = 1:length(disturbances)
    tic
x_only_pert = member_PO + disturbances(ii)/R*[1 0 0 0 0 0 ]';
x_only_pert(5) = find_ydot(member_PO,x_only_pert);
[~,X_xo,~,X_e_xo,~] = ode45(@Lagrange_CR3BP,[0 member_T*1000], x_only_pert, ...
    odeset('Events', @y_crossing_1k, 'RelTol', 3e-14, 'AbsTol', 1e-14), ...
    propatagor_opts);
% C = compute_C(x_only_pert(1:3),x_only_pert(4:6))

% plot(X_e_xo(:,2))
figure(x_only_p_map_x_z);
plot(X_e_xo(:,1),X_e_xo(:,3),'x','color', colors{ii})
axis equal
hold on
drawnow
figure(x_only_p_map_x_vx);
plot(X_e_xo(:,1),X_e_xo(:,4),'x','color', colors{ii})
hold on
drawnow
toc
% plot(X_e(:,1),X_e(:,6),'x')
end


figure(x_only_p_map_x_z);
title('Poincare Map, Example DRO x vs z')
x_label('x, DU'); y_label('z, DU')


figure(x_only_p_map_x_vx);
title('Poincare Map, Example DRO x vs v_x')
x_label('x, DU'); y_label('v_x, DU/TU')

figure
plot(X_xo(:,1),X_xo(:,2))
title('Example DRO With x-Only Perturbation')
x_label('x, DU'); y_label('y, DU')


%% Perturb Z
Z_only_p_map = figure(Position', hw_pub.figPosn);
z_only_p_map_z_vz = figure(Position', hw_pub.figPosn);
for ii = 1:length(disturbances)
    tic
Z_only_pert = member_PO + disturbances(ii)/R*[0 0 1 0 0 0 ]';
Z_only_pert(5) = find_ydot(member_PO,Z_only_pert);
[~,X_zo,~,X_e_zo,~] = ode45(@Lagrange_CR3BP,[0 member_T*1000], Z_only_pert, ...
    odeset('Events', @y_crossing_1k, 'RelTol', 3e-14, 'AbsTol', 1e-14), ...
    propatagor_opts);

figure(Z_only_p_map);
plot(X_e_zo(:,1),X_e_zo(:,3),'x','color', colors{ii})
axis equal; hold on; drawnow;
figure(z_only_p_map_z_vz);
plot(X_e_zo(:,3),X_e_zo(:,6),'x','color', colors{ii})
axis equal; hold on; drawnow;
toc
if disturbances(ii) == 1000
    figure(Position', hw_pub.figPosn);
    plot3(X_zo(:,1),X_zo(:,2),X_zo(:,3))
end
end
%% Perturb VX
vx_only_pert = member_PO + 0.075/R*P/(2*pi)*[0 0 0 1 0 0 ]';
vx_only_pert(5) = find_ydot(member_PO,vx_only_pert);
[~,~,~,X_e_vxo,~] = ode45(@Lagrange_CR3BP,[0 member_T*1000], vx_only_pert, ...
    odeset('Events', @y_crossing_1k, 'RelTol', 3e-14, 'AbsTol', 1e-14), ...
    propatagor_opts);

vx_only_p_map = figure(Position', hw_pub.figPosn);
plot(X_e_vxo(:,1),X_e_vxo(:,3),'x')
axis equal
vx_only_p_map_x_vx = figure(Position', hw_pub.figPosn);
plot(X_e_vxo(:,1),X_e_vxo(:,4),'x')
axis equal
% plot(X_e(:,1),X_e(:,6),'x')
%% Perturb XZ
xz_pert = member_PO + 1000/R/sqrt(2)*[1 0 1 0 0 0 ]';
xz_pert(5) = find_ydot(member_PO,xz_pert);
[~,~,~,X_e_xz,~] = ode45(@Lagrange_CR3BP,[0 member_T*1000], xz_pert, ...
    odeset('Events', @y_crossing_1k, 'RelTol', 3e-14, 'AbsTol', 1e-14), ...
    propatagor_opts);

xz_p_map = figure(Position', hw_pub.figPosn);
% plot(X_e_xo(:,2))
plot(X_e_xz(:,1),X_e_xz(:,3),'x')
axis equal
xz_p_map_x_vx = figure(Position', hw_pub.figPosn);
plot(X_e_xz(:,1),X_e_xz(:,4),'x')
axis equal