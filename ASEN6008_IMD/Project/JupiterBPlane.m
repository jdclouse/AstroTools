%% B Plane for all gravity assists

fprintf(['Launch: ' getDate(Launch_date) '\n'])
fprintf(['VGA: ' getDate(VGA_date) '\n'])
fprintf(['EGA1: ' getDate(EGA1_date) '\n'])
fprintf(['EGA2: ' getDate(EGA2_date) '\n'])
fprintf(['JOI: ' getDate(JGA_date) '\n'])
fprintf('\n')
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

fprintf('C3: '); disp(c3_store(use_traj(1), VGA_date_idx)); fprintf('\b\b km^2/s^2\n');
fprintf('V_inf JOI: '); disp(dv_final_store(EGA_date_idx, JGA_date_idx)); fprintf('\b\b km/s\n');

%% The target launch parameters
launch_C3 = lambert_out(1).lw_c3_store(...
    (Launch_date_idx),(VGA_date_idx));
launch_DLA = asin(Launch_v_inf_out(3)/norm(Launch_v_inf_out));
launch_RLA = atan2(Launch_v_inf_out(2),Launch_v_inf_out(1));
fprintf('Launch Targets:\n')
fprintf('C3 = %.6f km^2/s^2\n',launch_C3);
fprintf('DLA = %.6f deg\n',launch_DLA*180/pi);
fprintf('RLA = %.6f deg\n',launch_RLA*180/pi);
