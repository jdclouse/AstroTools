function [ value, isterminal, direction ] = detect_rp( t,X,r_lim )
%UNTITLED Stop when the SC apoapsis crosses target sma

CelestialConstants;

r = Euler2DCM('3',X(2))*[X(1);0;0];
v = Euler2DCM('3',X(2))*[X(3); ...
                        (X(1)*X(4));...
                         0];
[a,e,~,~,~,~] = cart2OE(r,v,Sun.mu);
rp = a*(1-e);
value = r_lim - rp;% Polar
isterminal = 1;
direction = 0;
end