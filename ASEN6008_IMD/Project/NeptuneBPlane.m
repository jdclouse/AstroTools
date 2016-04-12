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


%%

SGA_date_idx = use_traj(5);
NOI_date_idx = use_traj(6);
SGA_date = SGA_arr(SGA_date_idx);
NOI_date = NOI_arr(NOI_date_idx);
[r_saturn_SGA, v_saturn_SGA] = ...
    MeeusEphemeris(Saturn, SGA_date,Sun);
[r_neptune_NOI, v_neptune_NOI] = ...
    MeeusEphemeris(Neptune, NOI_date,Sun);

% The outgoing velocity on JGA.
[JGA_v_helio_out, SGA_v_helio_in] = lambert( r_jupiter_JGA, r_saturn_SGA, ...
    (SGA_date-JGA_date)*day2sec, ...
    1, Sun);
JGA_v_inf_out = JGA_v_helio_out - v_jupiter_JGA;
SGA_v_inf_in = SGA_v_helio_in - v_saturn_SGA;

% The outgoing velocity on SGA.
[SGA_v_helio_out, NOI_v_helio_in] = lambert( r_saturn_SGA, r_neptune_NOI, ...
    (NOI_date-SGA_date)*day2sec, ...
    1, Sun);
SGA_v_inf_out = SGA_v_helio_out - v_saturn_SGA;
NOI_v_inf_in = NOI_v_helio_in - v_neptune_NOI;


% JGA
[b_JGA, B_hat_JGA, B_plane_JGA, psi_JGA, rp_JGA] = ...
    BPlaneTarget(Jup_v_inf_in, JGA_v_inf_out, Jupiter.mu);
BT_JGA = dot(b_JGA*B_hat_JGA, B_plane_JGA(:,2));
BR_JGA = dot(b_JGA*B_hat_JGA, B_plane_JGA(:,3));
fprintf('JGA r_p = '); disp(rp_JGA);fprintf('\b\b km\n')
fprintf('        = '); disp(rp_JGA/Jupiter.R);fprintf('\b\b planetary radii\n')
fprintf('JGA turning angle = '); disp(psi_JGA*180/pi);fprintf('\b\b deg\n')
fprintf('JGA BT = '); disp(BT_JGA);fprintf('\b\b km\n')
fprintf('JGA BR = '); disp(BR_JGA);fprintf('\b\b km\n\n')

% SGA
[b_SGA, B_hat_SGA, B_plane_SGA, psi_SGA, rp_SGA] = ...
    BPlaneTarget(SGA_v_inf_in, SGA_v_inf_out, Saturn.mu);
BT_SGA = dot(b_SGA*B_hat_SGA, B_plane_SGA(:,2));
BR_SGA = dot(b_SGA*B_hat_SGA, B_plane_SGA(:,3));
fprintf('SGA r_p = '); disp(rp_SGA);fprintf('\b\b km\n')
fprintf('        = '); disp(rp_SGA/Saturn.R);fprintf('\b\b planetary radii\n')
fprintf('SGA turning angle = '); disp(psi_SGA*180/pi);fprintf('\b\b deg\n')
fprintf('SGA BT = '); disp(BT_SGA);fprintf('\b\b km\n')
fprintf('SGA BR = '); disp(BR_SGA);fprintf('\b\b km\n\n')