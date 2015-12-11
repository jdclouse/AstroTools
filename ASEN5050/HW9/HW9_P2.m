%% HW9 Problem 2
fprintf('\n');
clearvars -except function_list hw_pub toolsPath 
close all
CelestialConstants; % import useful constants

rp = 200+Mars.R;

% Transfer properties
a_xfer = (Earth.a + Mars.a)/2;
v_xfer_f = sqrt(2*Sun.mu/Mars.a - Sun.mu/a_xfer);

v_mars = sqrt(Sun.mu/Mars.a);
v_inf_mars = abs(v_mars-v_xfer_f);

v_assist_mars_max = abs(v_inf_mars + v_mars);
v_assist_mars_min = abs(v_inf_mars - v_mars);

fprintf(['Minimum heliocentric velocity from Mars gravity assist: %.3f'...
    ' km/s\n'],v_assist_mars_min)
fprintf(['Maximum heliocentric velocity from Mars gravity assist: %.3f'...
    ' km/s\n'],v_assist_mars_max)