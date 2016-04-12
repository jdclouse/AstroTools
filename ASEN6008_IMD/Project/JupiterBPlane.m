%% B Plane for all gravity assists
fprintf('B-plane targeting:\n')
% VGA
[b_VGA, B_hat_VGA, B_plane_VGA, psi_VGA, rp_VGA] = ...
    BPlaneTarget(VGA_v_inf_in, VGA_v_inf_out, Venus.mu);
BT_VGA = dot(b_VGA*B_hat_VGA, B_plane_VGA(:,2));
BR_VGA = dot(b_VGA*B_hat_VGA, B_plane_VGA(:,3));
fprintf('VGA r_p = '); disp(rp_VGA);fprintf('\b\b km\n')
fprintf('VGA turning angle = '); disp(psi_VGA*180/pi);fprintf('\b\b deg\n')
fprintf('VGA BT = '); disp(BT_VGA);fprintf('\b\b km\n')
fprintf('VGA BR = '); disp(BR_VGA);fprintf('\b\b km\n\n')

% EGA1
[b_EGA1, B_hat_EGA1, B_plane_EGA1, psi_EGA1, rp_EGA1] = ...
    BPlaneTarget(EGA1_v_inf_in, EGA1_v_inf_out, Earth.mu);
BT_EGA1 = dot(b_EGA1*B_hat_EGA1, B_plane_EGA1(:,2));
BR_EGA1 = dot(b_EGA1*B_hat_EGA1, B_plane_EGA1(:,3));
fprintf('EGA1 r_p = '); disp(rp_EGA1);fprintf('\b\b km\n')
fprintf('EGA1 turning angle = '); disp(psi_EGA1*180/pi);fprintf('\b\b deg\n')
fprintf('EGA1 BT = '); disp(BT_EGA1);fprintf('\b\b km\n')
fprintf('EGA1 BR = '); disp(BR_EGA1);fprintf('\b\b km\n\n')

% EGA2
[b_EGA2, B_hat_EGA2, B_plane_EGA2, psi_EGA2, rp_EGA2] = ...
    BPlaneTarget(EGA2_v_inf_in, EGA2_v_inf_out, Earth.mu);
BT_EGA2 = dot(b_EGA2*B_hat_EGA2, B_plane_EGA2(:,2));
BR_EGA2 = dot(b_EGA2*B_hat_EGA2, B_plane_EGA2(:,3));
fprintf('EGA2 r_p = '); disp(rp_EGA2);fprintf('\b\b km\n')
fprintf('EGA2 turning angle = '); disp(psi_EGA2*180/pi);fprintf('\b\b deg\n')
fprintf('EGA2 BT = '); disp(BT_EGA2);fprintf('\b\b km\n')
fprintf('EGA2 BR = '); disp(BR_EGA2);fprintf('\b\b km\n\n')
