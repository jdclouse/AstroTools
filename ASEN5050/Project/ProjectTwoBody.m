function [X_dot, sail_accel, fpa, phase_angle] = ProjectTwoBody(t,X,ctrl_law)


CelestialConstants; % import useful constants

W = 1361; %W/m2 = kg*m2/s3/m2
W = 0; %W/m2 = kg*m2/s3/m2
m_s = 40;
m_p = 500; %kg
m = m_p+m_s;
L_boom = 30;
h = L_boom*sin(pi/4)*2;
r = 0.88;
s = 0.94;
Bf = 0.79;
Bb = 0.55;
ef = 0.05;
eb=0.55;
rho_s = r*s;
rho_d = (Bf*r*(1-s)+(ef*Bf-eb*Bb)/(ef+eb))*3/2;
A_sail = h*h;

X = real(X);

beta = 0.05;

[~,e,~,~,~,f] = cart2OE(X(1:3),X(4:6),Sun.mu);

fpa = pi/2 - acos(dot(X(1:3),X(4:6))/(norm(X(1:3))*norm(X(4:6))));
fpa = (e*sin(f))/(1+e*cos(f));
    
phase_angle = atan2(X(2),X(1));

SRP_accel_STW = @(X) [cos(X(7))*cos(X(7))*cos(X(7));...
cos(X(7))*cos(X(7))*sin(X(7));0]*beta*Sun.mu/(norm(X(1:3))*norm(X(1:3)));

% Anonymous function to calculate SRP
SRP = @(X) W/speed_of_light/(norm(X(1:3)/au2km));

% Calculate body acceleration
SRP_accel_body = @(X) SRP(X)*A_sail/m*[...
    (1+rho_s)*cos(X(7))*cos(X(7))+2/3*rho_d*cos(X(7));...
    (1-rho_s)*cos(X(7))*sin(X(7)); 0]; % return 3-vec


% Acceleration in inertial frame
SRP_accel_inrtl = @(X,t) ...
    Euler2DCM('3',phase_angle)...
    *SRP_accel_STW(X); % return 3-vec

grav_accel = [-Sun.mu*X(1)/norm(X(1:3))^3;...
    -Sun.mu*X(2)/norm(X(1:3))^3;...
    -Sun.mu*X(3)/norm(X(1:3))^3];

X_dot = [X(4);X(5);X(6);...
    grav_accel; 0];

if strcmp(ctrl_law,'Yes')
    sail_accel = SRP_accel_inrtl(X,t);
    sail_accel = 1e-7*X(1:3)/norm(X(1:3));
elseif strcmp(ctrl_law,'Optimal')
    
    a_calc = pi/2-fpa;
    if imag(fpa)
    fpa
    t
X
imag(X)
norm(X(1:3))
norm(X(4:6))
a_calc = real(a_calc);
a_calc = max(a_calc, pi/2);
end
    X(7) = atan2(-3+sqrt(9+8*tan(a_calc)^2),4*tan(a_calc));
%     max_opt_sun_angle = atan(1/sqrt(2));
%     if X(7) > max_opt_sun_angle
%         X(7) = max_opt_sun_angle;
%     elseif X(7) < -max_opt_sun_angle
%         X(7) = -max_opt_sun_angle;
%     end
    sail_accel = SRP_accel_inrtl(X,t);
else
    sail_accel = zeros(3,1);
end

X_dot = X_dot + [0;0;0;sail_accel;0]; % return state-sized vec
if (imag(X_dot))
fprintf('Well theres your problem.')
disp X_dot
end