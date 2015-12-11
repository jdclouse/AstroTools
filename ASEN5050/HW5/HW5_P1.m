%% HW5 Problem 1
fprintf('\n');
clearvars -except function_list hw_pub toolsPath 
close all
CelestialConstants; % import useful constants

% quick function to compute velocity on the fly:
visviva = @(r,a) sqrt(2*Sun.mu/(r) - Sun.mu/a);

%% Earth->Venus
fprintf('Earth->Venus\n');
earth_v = visviva(Earth.a, Earth.a);
Venus_v = visviva(Venus.a, Venus.a);
a_xfer = (Earth.a + Venus.a)/2;
dv1 = visviva(Earth.a, a_xfer) - earth_v
dv2 = Venus_v - visviva(Venus.a, a_xfer)
dv_tot = abs(dv1) + abs(dv2)
T = pi*sqrt(a_xfer^3/Sun.mu)/day2sec/365.25

%% Earth->Mars
fprintf('Earth->Mars\n');
earth_v = visviva(Earth.a, Earth.a);
Mars_v = visviva(Mars.a, Mars.a);
a_xfer = (Earth.a + Mars.a)/2;
dv1 = visviva(Earth.a, a_xfer) - earth_v
dv2 = Mars_v - visviva(Mars.a, a_xfer)
dv_tot = abs(dv1) + abs(dv2)
T = pi*sqrt(a_xfer^3/Sun.mu)/day2sec/365.25

%% Earth->Jupiter
fprintf('Earth->Jupiter\n');
earth_v = visviva(Earth.a, Earth.a);
Jupiter_v = visviva(Jupiter.a, Jupiter.a);
a_xfer = (Earth.a + Jupiter.a)/2;
dv1 = visviva(Earth.a, a_xfer) - earth_v
dv2 = Jupiter_v - visviva(Jupiter.a, a_xfer)
dv_tot = abs(dv1) + abs(dv2)
T = pi*sqrt(a_xfer^3/Sun.mu)/day2sec/365.25