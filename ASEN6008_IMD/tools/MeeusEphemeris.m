function [ r, v ] = MeeusEphemeris( planet, JD , Sun)
%MeeusEphemeris Calculate planetary ephemeris. Works with
%CelestialConstants.m file
%   Outputs PV in km, km/s
% fcnPrintQueue(mfilename('fullpath')) % Add this code to code app 

T = (JD - 2451545)/36525;

if length(planet.Meeus.J200.a) == 1
    a = planet.Meeus.J200.a;%*au2km;
else
    T_pow = 1;
    a = 0;
    for ii = 1:length(planet.Meeus.J200.a) 
        a = a + planet.Meeus.J200.a(ii)*T_pow;
        T_pow = T_pow*T;
    end
end

L = planet.Meeus.J200.L(1) ...
    + planet.Meeus.J200.L(2)*T ...
    + planet.Meeus.J200.L(3)*T*T ...
    + planet.Meeus.J200.L(4)*T*T*T;
e = planet.Meeus.J200.e(1) ...
    + planet.Meeus.J200.e(2)*T ...
    + planet.Meeus.J200.e(3)*T*T ...
    + planet.Meeus.J200.e(4)*T*T*T;

i = planet.Meeus.J200.i(1) ...
    + planet.Meeus.J200.i(2)*T ...
    + planet.Meeus.J200.i(3)*T*T ...
    + planet.Meeus.J200.i(4)*T*T*T;

RAAN = planet.Meeus.J200.RAAN(1) ...
    + planet.Meeus.J200.RAAN(2)*T ...
    + planet.Meeus.J200.RAAN(3)*T*T ...
    + planet.Meeus.J200.RAAN(4)*T*T*T;

Pi = planet.Meeus.J200.Pi(1) ...
    + planet.Meeus.J200.Pi(2)*T ...
    + planet.Meeus.J200.Pi(3)*T*T ...
    + planet.Meeus.J200.Pi(4)*T*T*T;

% Convert everything to radians!
L = L*pi/180;
i = i*pi/180;
RAAN = RAAN*pi/180;
Pi = Pi*pi/180;

w = Pi - RAAN;

M = L - Pi;

e2 = e*e;
e3 = e*e2;
e4 = e*e3;
e5 = e*e4;

C_cen = (2*e-e3/4+5/96*e5)*sin(M) + (5/4*e2-11/24*e4)*sin(2*M) ...
    + (13/12*e3-43/64*e5)*sin(3*M) + 103/96*e4*sin(4*M) ...
    + 1097/960*e5*sin(5/M);

f = M + C_cen;

[r, v] = OE2cart(a, e, i, RAAN, w, f, Sun.mu);

end

