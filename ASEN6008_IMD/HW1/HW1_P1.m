%% IMD HW 1, Problem 1
% John Clouse
%% Initialize
clearvars -except hw_pub function_list

h_i = 400; %km
h_f = 400; %km

mu_sun = 1.327e11; %km3/s2
mu_earth = 3.986e5; %km3/s2
mu_mars = 4.305e4; %km3/s2

km_per_AU = 1.4959787e8; %km
a_earth = 1*km_per_AU; %km
a_mars = 1.52367934*km_per_AU; %km 
r_earth = 6378.1363; %km
r_mars = 3397.2; %km

%% Calculate transfer orbit, velocities
a_xfer = (a_earth + a_mars)/2;
v_earth = sqrt(mu_sun/a_earth);
v_mars = sqrt(mu_sun/a_mars);

v_xfer_i = sqrt(2*mu_sun/a_earth - mu_sun/a_xfer);
v_xfer_f = sqrt(2*mu_sun/a_mars - mu_sun/a_xfer);
fprintf('Departure velocity: %.3f km/s\n',v_xfer_i);
fprintf('Arrival velocity: %.3f km/s\n',v_xfer_f);

v_esc_earth = sqrt(2*mu_earth/(r_earth + h_i));
v_park_earth = sqrt(mu_earth/(r_earth + h_i));

v_esc_mars = sqrt(2*mu_mars/(r_mars + h_f));
v_park_mars = sqrt(mu_mars/(r_mars + h_f));

v_inf_dep = v_xfer_i-v_earth;
v_inf_arr = v_mars-v_xfer_f;
v_hyp_dep = sqrt(2*mu_earth/(r_earth + h_i) + v_inf_dep^2);
v_hyp_arr = sqrt(2*mu_mars/(r_mars + h_f) + v_inf_arr^2);
dv_departure = v_hyp_dep - v_park_earth;
dv_arrival = v_hyp_arr - v_park_mars;


fprintf('Departure dV: %.3f km/s\n',dv_departure);
fprintf('Arrival dV: %.3f km/s\n',dv_arrival);
fprintf('\t Initial velocity = %.3f km/s\n', v_park_earth);
fprintf('\t Hyperbolic departure velocity = %.3f km/s\n', v_hyp_dep);
fprintf('\t Departure dV = %.3f km/s\n', v_hyp_dep-v_park_earth);
fprintf('\t Hyperbolic arrival velocity = %.3f km/s\n', v_hyp_arr);
fprintf('\t Final Mars velocity = %.3f km/s\n', v_park_mars);
fprintf('\t Arrival dV = %.3f km/s\n', v_hyp_arr-v_park_mars);
T = pi*sqrt(a_xfer^3/mu_sun)/3600/24;
fprintf('Transfer time: %.3f days\n', T);