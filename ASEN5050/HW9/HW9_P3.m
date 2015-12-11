%% HW9 Problem 3
fprintf('\n');
clearvars -except function_list hw_pub toolsPath 
close all
CelestialConstants; % import useful constants

rp = Jupiter.R*2;

% Transfer properties
a_xfer = (Earth.a + Jupiter.a)/2;
v_xfer_f = sqrt(2*Sun.mu/Jupiter.a - Sun.mu/a_xfer);

v_jup = sqrt(Sun.mu/Jupiter.a);
v_inf_jup = abs(v_jup-v_xfer_f);

v_assist_jup_max = abs(v_inf_jup + v_jup);

v_esc_jup = sqrt(2*Sun.mu/Jupiter.a);

spec_energy = v_assist_jup_max^2/2-Sun.mu/Jupiter.a; %J/kg

fprintf('Max velocity from gravity assist: %.3f km/s\n', v_assist_jup_max)
fprintf('Heliocentric escape velocity at Jupiter: %.3f km/s\n', v_esc_jup)
fprintf('The spacecraft can escape the solar system in this mission plan.\n')
fprintf('Max specific energy: %.3f J/kg\n', spec_energy)