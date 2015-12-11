%% Problem 1
clear
a = 7078.1363*(1-0.22)
n = sqrt(398600.4415/a/a/a)
E1 = atan2(sind(50)*sqrt(1-0.22*0.22),0.22+cosd(50))
E2 = atan2(sind(230)*sqrt(1-0.22*0.22),0.22+cosd(230)) + 2*pi
M1 = E1 - 0.22*sin(E1)
M2 = E2 - 0.22*sin(E2)
delta_time = (M2-M1)/n
delta_time = delta_time/3600

%% Problem 2
clear
a = 8580;
e = 0.39;
i = 150.9;
raan = 275.0;
w = 110.1;
f1 = 230;
rp = 8580*(1-0.39)
p = 8580*(1-0.39*0.39)
f_imp = acosd((p/6378.1363-1)/e)
f_imp = 360-f_imp
n = sqrt(398600.4415/a^3)
E1 = atan2(sind(f1)*sqrt(1-e*e),e+cosd(f1)) + 2*pi
E_imp = atan2(sind(f_imp)*sqrt(1-e*e),e+cosd(f_imp)) + 2*pi
M1 = E1 - e*sin(E1)
M_imp = E_imp - e*sin(E_imp)
delta_time = (M_imp - M1)/n
delta_time = delta_time/3600
v_imp = sqrt(2*398600.4415/6378.1363 - 398600.4415/a)

%% Problem 3
clear
% a)
w = sqrt(398600.4415/(6378.1363+380)^3)
t2 = 60*8
x_b2 = 2*0.2/w*(1 - cos(w*t2))
y_b2 = 0.2*(4/w*sin(w*t2) - 3*t2)
z_a2 = 0.4/w*sin(w*t2)
% b)
t3 = 60*20
x_b3 = 2*0.2/w*(1 - cos(w*t3))
y_b3 = 0.2*(4/w*sin(w*t3) - 3*t3)
t_elapsed = 12*60
A = [sin(w*t_elapsed)/w  2/w*(1-cos(w*t_elapsed)), 0;...
    2/w*cos(w*t_elapsed)-2/w, 4/w*sin(w*t_elapsed)-3*t_elapsed, 0;...
    0, 0, sin(w*t_elapsed)/w];
v_imp_1_plus = inv(A)*[x_b3; y_b3; -z_a2*cos(w*t_elapsed)]
zdot_a2 = 0.4*cos(w*t2)
v_burnA = v_imp_1_plus-[0;0;zdot_a2]
vb_2 = [2*.2*sin(w*t2); 0.2*(4*cos(w*t2)-3);0]
v_burnA_wrtB = v_burnA-vb_2
vb_3 = [2*.2*sin(w*t3); 0.2*(4*cos(w*t3)-3);0]
% c)
va_prerdv = [v_imp_1_plus(1)*cos(w*t_elapsed) + ...
    2*v_imp_1_plus(2)*sin(w*t_elapsed);...
    -2*v_imp_1_plus(1)*sin(w*t_elapsed) + ...
    v_imp_1_plus(2)*(4*cos(w*t_elapsed)-3);...
    -z_a2*sin(w*t_elapsed) + v_imp_1_plus(3)*cos(w*t_elapsed)]

%% Problem 4
clear
theta_dot_earth = 360*.99726968/24/3600
P = 25/theta_dot_earth
a = ((P/2/pi)^2*398600.4415)^(1/3)
h = a - 6378.1363

%% Problem 5
clear 
% a)
a = (362600 + 405400)/2
P = 2*pi*sqrt(a^3/398600.4415)
P = P/3600/24
% b)
e = (405400 - 362600)/(362600 + 405400)
p = a*(1-e*e)
f = acosd((p/380000-1)/e)
n = sqrt(398600.4415/a^3)
E = atan2(sind(f)*sqrt(1-e*e),e+cosd(f))
M = E - e*sin(E)
delta_time = 2*M/n
delta_time = delta_time/3600/24
percentage = delta_time/P*100
%% Problem 6
clear
r = [8800; -1100; 5500];
norm(r)
phi = asind(r(3)/norm(r))
theta_sat = atan2d(r(2),r(1))
theta_LST = 90 - 105.27
R_LST = [cosd(theta_LST) sind(theta_LST) 0;
    -sind(theta_LST) cosd(theta_LST) 0;
    0 0 1];
R_lat = [cosd(90-40.015) 0 -sind(90-40.015);
    0 1 0;
    sind(90-40.015) 0 cosd(90-40.015)];
r_SEZ = R_lat*R_LST*r
range = r_SEZ - [0;0;6379.77]
el = asind(range(3)/norm(range))
az = 180 - atan2d(range(2),range(1))

%% Problem 7
clear
2*asind((sqrt(3)-1/sqrt(3)-1)/2)