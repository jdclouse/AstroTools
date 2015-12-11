%% HW9 Problem 1
fprintf('\n');
clearvars -except function_list hw_pub toolsPath 
close all
CelestialConstants; % import useful constants

h_earth = 185; %km
h_mars = 300; %km

r_soi_earth = Earth.a*(Earth.m/Sun.m)^(2/5)
r_soi_mars = Mars.a*(Mars.m/Sun.m)^(2/5)
r_soi_Earth_Moon = Earth.a*((Moon.m+Earth.m)/Sun.m)^(2/5)

% Transfer properties
a_xfer = (Earth.a + Mars.a)/2;
v_xfer_i = sqrt(2*Sun.mu/Earth.a - Sun.mu/a_xfer);
v_xfer_f = sqrt(2*Sun.mu/Mars.a - Sun.mu/a_xfer);

v_earth = sqrt(Sun.mu/Earth.a);
v_mars = sqrt(Sun.mu/Mars.a);

v_inf_earth = abs(v_xfer_i-v_earth);
v_inf_mars = abs(v_mars-v_xfer_f);

v_park_earth = sqrt(Earth.mu/(Earth.R+h_earth));
v_park_mars = sqrt(Mars.mu/(Mars.R+h_mars));

dv_earth_inj = sqrt(2*Earth.mu/(Earth.R+h_earth)+v_inf_earth^2)...
    -v_park_earth;
dv_mars_ins = sqrt(2*Mars.mu/(Mars.R+h_mars)+v_inf_mars^2)-v_park_mars;

T_xfer = 2*pi*sqrt(a_xfer^3/Sun.mu)/2/3600/24;


fprintf('a) Earth SOI: %.0f km\n',r_soi_earth)
fprintf('   Mars SOI: %.0f km\n',r_soi_mars)
fprintf('   Earth-Moon SOI: %.0f km\n',r_soi_Earth_Moon)
fprintf('b) Earth heliocentric departure velocity: %.3f km/s\n',v_xfer_i)
fprintf('   Mars heliocentric arrival velocity: %.3f km/s\n',v_xfer_f)
fprintf('c) Earth departure v_inf: %.3f km/s\n',v_inf_earth)
fprintf('   Mars arrival v_inf: %.3f km/s\n',v_inf_mars)
fprintf('d) Mars transfer injection dv: %.3f km/s\n',dv_earth_inj)
fprintf('e) Mars insertion dv: %.3f km/s\n',dv_mars_ins)
fprintf('f) Transfer time: %.1f days\n',T_xfer)