%% Resonant Orbit
% Constructing a 2:1 resonant orbit
%% Ephemerides and v_infs
good_reso = 0;

use_traj = best_traj;
% use_traj = lowest_vf_traj;
% use_traj = lowest_C3_traj;
Launch_date = Launch_dep(use_traj(1));
Launch_date_idx = use_traj(1);
VGA_date_idx = use_traj(2);
EGA_date_idx = use_traj(3);
JGA_date_idx = use_traj(4);
VGA_date = VGA_arr(VGA_date_idx);
EGA1_date = EGA1_arr(EGA_date_idx);
EGA2_date = EGA2_arr(EGA_date_idx);
JGA_date = JGA_arr(JGA_date_idx);

incoming_v = eval(inc_vel(Venus_Earth));
outgoing_v = eval(out_vel(Earth_Jupiter));
% EGA1_v_inf = lambert_out(2).long_way_dv2_store(VGA_date_idx, EGA_date_idx);
% EGA2_v_inf = lambert_out(3).long_way_dv1_store( EGA_date_idx, JGA_date_idx);
EGA1_v_inf = incoming_v(VGA_date_idx, EGA_date_idx);
EGA2_v_inf = outgoing_v(EGA_date_idx, JGA_date_idx);

[r_earth_launch, v_earth_launch] = ...
    MeeusEphemeris(Earth, Launch_date,Sun);
[r_venus_vga, v_venus_vga] = ...
    MeeusEphemeris(Venus, VGA_date,Sun);
[r_earth_ega1, v_earth_ega1] = ...
    MeeusEphemeris(Earth, EGA1_date,Sun);
[r_earth_ega2, v_earth_ega2] = ...
    MeeusEphemeris(Earth, EGA2_date,Sun);
[r_jupiter_JGA, v_jupiter_JGA] = ...
    MeeusEphemeris(Jupiter, JGA_date,Sun);

Launch_to_VGA_direction = -1;

% Launch to VGA
[Launch_v_helio_out, VGA_v_helio_in] = lambert( r_earth_launch, r_venus_vga, ...
    (VGA_date-Launch_date)*day2sec, ...
    Earth_Venus.lambert, Sun, 0, 1e-6);
Launch_v_inf_out = Launch_v_helio_out - v_earth_launch;
VGA_v_inf_in = VGA_v_helio_in - v_venus_vga;

% Incoming velocity on EGA1
[VGA_v_helio_out, EGA1_v_helio_in] = lambert( r_venus_vga, r_earth_ega1, ...
    (EGA1_date-VGA_date)*day2sec, ...
    Venus_Earth.lambert, Sun);
VGA_v_inf_out = VGA_v_helio_out - v_venus_vga;
EGA1_v_inf_in = EGA1_v_helio_in - v_earth_ega1;

% The outgoing velocity on EGA2.
[EGA2_v_helio_out, JOI_v_helio] = lambert( r_earth_ega2, r_jupiter_JGA, ...
    (JGA_date-EGA2_date)*day2sec, ...
    Earth_Jupiter.lambert, Sun, 0, 1e-6);
EGA2_v_inf_out = EGA2_v_helio_out - v_earth_ega2;
Jup_v_inf_in = JOI_v_helio - v_jupiter_JGA;

%%
[b_VGA, B_hat_VGA, B_plane_VGA, psi_VGA, rp_VGA] = ...
    BPlaneTarget(VGA_v_inf_in, VGA_v_inf_out, Venus.mu);
BT_VGA = dot(b_VGA*B_hat_VGA, B_plane_VGA(:,2));
BR_VGA = dot(b_VGA*B_hat_VGA, B_plane_VGA(:,3));
% fprintf('VGA r_p = '); disp(rp_VGA);fprintf('\b\b km\n')
% fprintf('VGA turning angle = '); disp(psi_VGA*180/pi);fprintf('\b\b deg\n')
% fprintf('VGA BT = '); disp(BT_VGA);fprintf('\b\b km\n')
% fprintf('VGA BR = '); disp(BR_VGA);fprintf('\b\b km\n\n')
if rp_VGA < 6000
    fprintf('no good VGA\n');
%     rp_VGA
    return
end

%% The phi calculation
% period is 2 years
P = 2*365.242189*day2sec;
a_reso = (P*P/4/pi/pi*Sun.mu)^(1/3); %sma

% The velocity immediately after the first gravity assist
V_sc_sun = sqrt(Sun.mu*(2/norm(r_earth_ega1) - 1/a_reso));

theta = acos((-norm(V_sc_sun)^2 + norm(v_earth_ega1)^2 + EGA1_v_inf^2)...
    /(2*EGA1_v_inf*norm(v_earth_ega1)));

% The velocity is calculated in the Earth VNC frame. Calculate the
% transformation from VNC to ecliptic coords.
V_hat = v_earth_ega1/norm(v_earth_ega1);
N_hat = cross(r_earth_ega1, v_earth_ega1)...
    /norm(cross(r_earth_ega1, v_earth_ega1)); % angular momentum
C_hat = cross(V_hat, N_hat);
T_VNC2Ecl = [V_hat N_hat C_hat];

% Cycle through the locus of possible orbits
r_min = 7000; % Minimum radius of earth approach
r_min = Earth.R + 300; % Minimum radius of earth approach
cos_term = cos(pi-theta);
sin_term = sin(pi-theta);
acceptable_phi = [];
acceptable_radii = [];
acceptable_EGA1_out = [];
acceptable_EGA2_in = [];
angles = 0:0.01:2*pi;
psi1_store = zeros(1,length(angles));
psi2_store = zeros(1,length(angles));
rp1_store = zeros(1,length(angles));
rp2_store = zeros(1,length(angles));
cnt = 0;
for phi = angles
    cnt = cnt+1;
    V_GA1_out = T_VNC2Ecl*EGA1_v_inf...
        *[cos_term; sin_term*cos(phi);-sin_term*sin(phi)];
    
    V_GA2_in = V_GA1_out + v_earth_ega1 - v_earth_ega2;
    
    % Determine if there is impact potential.
    % GA1:
    % Turning angle
    psi1 = acos(dot(EGA1_v_inf_in,V_GA1_out)...
        /norm(EGA1_v_inf_in)/norm(V_GA1_out));
    psi1_store(cnt) = psi1;
    % Closes approach to planet
    rp1 = Earth.mu/(norm(V_GA1_out))^2*(1/cos((pi-psi1)/2)-1);
    rp1_store(cnt) = rp1;
    % GA2:
    % Turning angle
    psi2 = acos(dot(V_GA2_in,EGA2_v_inf_out)...
        /norm(V_GA2_in)/norm(EGA2_v_inf_out));
    psi2_store(cnt) = psi2;
    % Closes approach to planet
    rp2 = Earth.mu/(norm(EGA2_v_inf_out))^2*(1/cos((pi-psi2)/2)-1);
    rp2_store(cnt) = rp2;
    
    if rp1 > r_min && rp2 > r_min
        acceptable_phi = [acceptable_phi phi]; %#ok<AGROW>
        acceptable_radii = [acceptable_radii [rp1;rp2]]; %#ok<AGROW>
        acceptable_EGA1_out = [acceptable_EGA1_out V_GA1_out]; %#ok<AGROW>
        acceptable_EGA2_in = [acceptable_EGA2_in V_GA2_in]; %#ok<AGROW>
    end
end
% figure('Position', hw_pub.figPosn);
% hold on
% plot(angles,rp1_store)
% plot(angles,rp2_store, 'r')
% plot(angles,Earth.R*ones(1,length(angles)), '--k')
% xlabel('\phi (radians)')
% xlim([0 2*pi])
% ylabel('r_p (km)')
% title('Closest Earth Approach for Resonant E-E Orbit')
% plot(angles,psi1_store)
% plot(angles,psi2_store)
% I'm choosing the phi such that both passes are as far away as possible to
% minimize drag and other near-earth affects.
if max(acceptable_phi) > 0
    good_reso = 1;
[~,xxx] = min(max(acceptable_radii,[],2));
[~,max_r_idx] = max(acceptable_radii(xxx,:));
% fprintf('Earth Resonant Orbit:\n')
% fprintf('phi = ');disp(acceptable_phi(max_r_idx)*180/pi);
% fprintf('\b\b deg\n');
% fprintf('EGA1 r_p = ');disp(acceptable_radii(1,max_r_idx));
% fprintf('\b\b km\n');
% fprintf('EGA2 r_p = ');disp(acceptable_radii(2,max_r_idx));
% fprintf('\b\b km\n');
% fprintf('EGA1 V_out = \n');disp(acceptable_EGA1_out(:,max_r_idx));
% fprintf('\b\b km/s\n');
% fprintf('EGA2 V_in = \n');disp(acceptable_EGA2_in(:,max_r_idx));
% fprintf('\b\b km/s\n');
% fprintf('\n');
% EGA1_v_inf_out = acceptable_EGA1_out(:,max_r_idx);
% EGA2_v_inf_in = acceptable_EGA2_in(:,max_r_idx);
else
    fprintf('No valid Resonant Orbit!!!\n');
end

