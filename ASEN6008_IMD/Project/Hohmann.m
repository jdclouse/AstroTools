%%Hohmann
v_earth_circ = sqrt(Sun.mu/Earth.a);
v_jup_circ = sqrt(Sun.mu/Jupiter.a);

a_xfer = (Earth.a + Jupiter.a)/2;
vp_xfer = sqrt(2*Sun.mu/Earth.a - Sun.mu/a_xfer);
va_xfer = sqrt(2*Sun.mu/Jupiter.a - Sun.mu/a_xfer);

T = pi*sqrt(a_xfer^3/Sun.mu)/3600/24/365.25;
C3 = (vp_xfer-v_earth_circ)^2
v_inf_arr = va_xfer - v_jup_circ


v_inf_arr = 5.59;
rp = 4000000;
dv_ins = sqrt(Jupiter.mu*2/rp + v_inf_arr^2) - sqrt(Jupiter.mu/rp)

r_JSOI = Jupiter.a*(Jupiter.m/Sun.m)^(2/5)
r_Ganymede = 1071.1e3;

R_G = 2631;

v_jup_peri = sqrt(Jupiter.mu*2/r_Ganymede + v_inf_arr^2)
v_gan_circ = sqrt(Jupiter.mu/r_Ganymede);

v_inf_gan = v_jup_peri - v_gan_circ;

phi_s = acos(1/(1+v_inf_gan^2*(R_G+200)/9885));
phi = pi-2*phi_s

r_a_final = 1e7;
a_final = (r_a_final+r_Ganymede)/2;
P_gan_days = 7.155
a_final = ((20*P_gan_days*24*3600/2/pi)^2*Jupiter.mu)^(1/3); %12:1 resonance
r_a_final = 2*a_final-r_Ganymede;
v_p_final = sqrt(2*Jupiter.mu/r_Ganymede - Jupiter.mu/a_final);
dv_ins = sqrt(Jupiter.mu*2/r_Ganymede + v_inf_arr^2) - v_p_final
2*pi*sqrt(a_final^3/Jupiter.mu)/3600/24

%% Target the B plane
e_hyp = ((Jupiter.mu*2/r_Ganymede + v_inf_arr^2)-Jupiter.mu/r_Ganymede)*r_Ganymede/Jupiter.mu;
a_hyp = -Jupiter.mu/v_inf_arr/v_inf_arr;
B_mag = -a_hyp*sqrt(e_hyp*e_hyp-1)