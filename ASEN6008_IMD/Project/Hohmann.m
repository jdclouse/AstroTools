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

%% Europa SOI
r_ESOI = 671.4e3*(3196/Jupiter.mu)^(2/5)

EurInj = norm([       -0.2021268970271654 
                    -0.09229850003525951
                     0.04043692774225515])
              
EurIns = norm([         -5.12522143513952 
                  1.364242052659392e-015
                 5.471179065352771e-016])
              
%% Ganymede
Period:         142.3265208830074 day  Rad. Apo:   1.4706999999997567e+007 km 14 Feb 2026 12:00:00.004 UTCG
Period:         79.47088010196345 day  Rad. Apo:   9.6141675683632549e+006 km 26 Apr 2026 15:55:03.999 UTCG
Period:         43.39294109872095 day  Rad. Apo:   6.1076384538077535e+006 km 14 Aug 2026 18:23:50.746 UTCG
Period:          28.5459066380272 day  Rad. Apo:   4.4061968318778453e+006 km 14 Oct 2026 07:08:39.955 UTCG

JTCM1=norm([     -0.003992651148234995 
                      0.2392092755335887 
                    -0.04211976786582281]) 
              
JTCM2=norm([                   -0.002402016575019953 
                    -0.01221331039729084 
                    -0.07119955383848155 ])
              
JTCM3=norm([                   -0.005394208346032947 
                   -0.009965304490565962 
                    -0.09728181861921918 ])
                
v1 = sqrt(2*Jupiter.mu/r_Ganymede - Jupiter.mu/1.4706999999997567e+007)
v2 = sqrt(2*Jupiter.mu/r_Ganymede - Jupiter.mu/9.6141675683632549e+006)
v3 = sqrt(2*Jupiter.mu/r_Ganymede - Jupiter.mu/6.1076384538077535e+006)
v4 = sqrt(2*Jupiter.mu/r_Ganymede - Jupiter.mu/4.4061968318778453e+006)

v1 - v2
v2 - v3
v3 - v4
