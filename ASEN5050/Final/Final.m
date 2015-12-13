%% Problem 1
clear all
CelestialConstants

ws = sqrt((Earth.mu+Moon.mu)/(384000-(384000*Moon.m/(Earth.m+Moon.m)))^3)

% T_tot = 5/3*pi/ws

h_phasing=1000;
a_xfer1 = Earth.R+(300+h_phasing)/2;
t_xfer1 = 2*pi*sqrt(a_xfer1^3/Earth.mu)/2

n_p = sqrt(Earth.mu/(Earth.R+h_phasing)^3)

a_xferf = (Earth.R+h_phasing+384000)/2
t_xferf = 2*pi*sqrt(a_xferf^3/Earth.mu)/2
n_xfer = sqrt(Earth.mu/a_xferf^3)

% L4_init_phase = atan((Moon.m/(Earth.m+Moon.m)-1/2))
L4_phase_angle_during_init_transfer = t_xfer1*ws
L4_phase_angle_during_final_transfer = t_xferf*ws

revs = 2;
phase_time = (pi/3-pi+L4_phase_angle_during_init_transfer+...
    L4_phase_angle_during_final_transfer -3*pi + revs*2*pi)/(n_p-ws)

dv1 = sqrt(2*Earth.mu/(Earth.R+300)-Earth.mu/a_xfer1)...
    -sqrt(Earth.mu/(Earth.R+300))
dv2 = sqrt(Earth.mu/(Earth.R+1000))...
    -sqrt(2*Earth.mu/(Earth.R+1000)-Earth.mu/a_xfer1)
dv3 = sqrt(2*Earth.mu/(Earth.R+1000)-Earth.mu/a_xferf)...
    -sqrt(Earth.mu/(Earth.R+1000))
dv4 = 384000*ws - sqrt(2*Earth.mu/384000-Earth.mu/a_xferf)

dv_incl_max = 2*sqrt(Earth.mu/(Earth.R+300))*sind(33.145396/2)
% v_L4 = sqrt()
%% Problem 2
r = 6378.1363+1.655064; %km

JD = computeJD(2015,12,8,13-7,31,41);
T_UT1 = ((JD)-2451545)/36525;

GMST = 67310.54841 ...
    + (876600*3600+8640184.812866) * T_UT1 ...
    + 0.093104 * T_UT1 * T_UT1 ...
    - 6.2e-6 * T_UT1 * T_UT1 * T_UT1;
while GMST > 86400
GMST = GMST-86400
end
LST = GMST/240-105

phi_gd = atand(tand(40)/(1-Earth.oblate_ecc^2));
r_eci = Euler2DCM('23', [-(90-phi_gd), -LST]*pi/180)*[0;0;r]
v = r_eci/r*0.002+cross([0;0;Earth.spin_rate],r_eci)
% v = cross([0;0;Earth.spin_rate],r_ecf)

[a,e,i,w,RAAN,f] = cart2OE(r_eci,v,Earth.mu)
ra = a*(1+e)
rp = a*(1-e)
h_max = ra-r

%% Problem 3
clear
CelestialConstants
a = 3*149597870.7;
m_ast = 1e15;
SOI = (m_ast/1.9891e30)^(2/5)*a

G = Earth.mu/Earth.m;%6.674e-20;
mu_ast = G*m_ast
P = 2*pi*sqrt(40^3/mu_ast)
P/3600
v_circ = sqrt(mu_ast/40)
a_new = 1/(2/40-(v_circ-.001)^2/mu_ast)
e_xfer = 40/a_new-1
n_xfer = sqrt(mu_ast/a_new^3)
Mf = pi+n_xfer*6*3600
Mf*180/pi
f = E2f(M2E(Mf,e_xfer),e_xfer)
f*180/pi
rf = a_new*(1-e_xfer^2)/(1+e_xfer*cos(f))
fpa = atan2(e_xfer*sin(f),1+e_xfer*cos(f))
fpa*180/pi
v_circ_f = sqrt(mu_ast/rf)
v_xfer_f = sqrt(2*mu_ast/rf-mu_ast/a_new)
dV = [0;v_circ_f]-[sin(fpa);cos(fpa)]*v_xfer_f
norm(dV)


%% Problem 4
clear
CelestialConstants

r_venus = [48965315.1 96179438.8 0.0]'; %km
v_venus = [-31.322263 15.730492 0.0]'; %km

v_approach = [-28.123456 8.654321 0.0]'; %km

SOI = (Venus.m/Sun.m)^(2/5)*norm(r_venus)
energy = norm(v_approach)^2/2-Sun.mu/(norm(r_venus))
V_inf_app = v_approach-v_venus
norm(V_inf_app)
energy_venus = norm(V_inf_app)^2/2-Venus.mu/SOI

%ccw 
V_inf_dep_ccw = Euler2DCM('3',[42*pi/180])*V_inf_app
V_dep_ccw = V_inf_dep_ccw+v_venus
energy = norm(V_dep_ccw)^2/2-Sun.mu/(norm(r_venus))

%cw 
V_inf_dep_cw = Euler2DCM('3',[-42*pi/180])*V_inf_app
V_dep_cw = V_inf_dep_cw+v_venus
energy = norm(V_dep_cw)^2/2-Sun.mu/(norm(r_venus))

rp = Venus.mu/norm(V_inf_app)*(1/cosd((180-42)/2)-1)

%% Problem 5
clear
CelestialConstants
h = 185; %km
v_inf = 7-(sqrt(2*Earth.mu/(h+Earth.R)) - sqrt(Earth.mu/(h+Earth.R)))

v_earth = sqrt(Sun.mu/au2km)

vp = v_inf+v_earth
a_xfer1 = 1/(2/au2km-norm(vp)^2/Sun.mu)
a_xfer1/au2km
ra1 = 2*a_xfer1-au2km
ra1/au2km

%b
va = v_earth-v_inf
a_xfer2 = 1/(2/au2km-norm(va)^2/Sun.mu)
a_xfer2/au2km
rp2 = 2*a_xfer2-au2km
rp2/au2km

%c 
i = asind(v_inf/v_earth)

%% Problem 6
clear
CelestialConstants
rp = 7e3;
ra = 8e3;
i = 100;
P = 60*100; %s
R = 6e3;
r_planet = 2*au2km;

a = (rp+ra)/2
mu = 2*pi*a^3/P^2

P_planet = 2*pi*sqrt(r_planet^3/Sun.mu)
RAAN_dot_ss = 2*pi/P_planet 
e = (ra-rp)/(ra+rp)

J2 = -RAAN_dot_ss*sqrt(a^3/mu)*2*a^2*(1-e*e)^2/(3*R*R*cosd(i))