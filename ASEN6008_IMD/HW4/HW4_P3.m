%% John Clouse IMD HW 4 Problem 3
% 
%% Initialize
clearvars -except hw_pub function_list

CelestialConstants;

JD_launch = 2447807.5;
JD_Venus = 2447932.5;
JD_Earth1 = 2448235.5;

%% Planet positions and velocities, V_inf in and out
[r_earth_L,~] = MeeusEphemeris(Earth, JD_launch,Sun);
[r_venus, v_venus] = MeeusEphemeris(Venus,JD_Venus,Sun);
[r_earth_1,~] = MeeusEphemeris(Earth, JD_Earth1,Sun);

[~,v_in] = lambert(r_earth_L, r_venus,(JD_Venus-JD_launch)*day2sec,1,Sun);
[v_out,~] = lambert(r_venus, r_earth_1,(JD_Earth1-JD_Venus)*day2sec,-1,Sun);

v_inf_in = v_in - v_venus;
v_inf_out = v_out - v_venus;

%% Flyby params
psi = acos(dot(v_inf_in,v_inf_out)/norm(v_inf_in)/norm(v_inf_out));
rp = Earth.mu/(norm(v_inf_in)^2)*(1/cos((pi-psi)/2)-1);
hp = rp - Venus.R; %km

%% Energy
energy_pre_flyby = norm(v_in)^2/2 - Sun.mu/norm(r_venus);
energy_post_flyby = norm(v_out)^2/2 - Sun.mu/norm(r_venus);
percent_change = (energy_post_flyby - energy_pre_flyby)...
    /abs(energy_pre_flyby)*100;

%% Results
% The v-infinities do not match exactly because the events aren't precisely
% targeted. True event JDs can be found to make them line up even better.
% The altitude of closes approach to Venus is shown below. The energy
% changed by 15%, making a trajectory with a higher aphelion.
fprintf('Closest approch altitude: %.3f km\n',hp);
fprintf('Heliocentric energy before: %.3f km^2/s^2\n',energy_pre_flyby);
fprintf('Heliocentric energy after: %.3f km^2/s^2\n',energy_post_flyby);