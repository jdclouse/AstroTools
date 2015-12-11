function [ value, isterminal, direction ] = detectDistLimit( t,X,dist )
%UNTITLED7 Stop when the radial vector crosses Mars Orbit
%   Detailed explanation goes here
CelestialConstants;

% value = Mars.a - norm(real(X(1:3))); %Cartesian
value = dist - X(1);% Polar
isterminal = 1;
direction = 0;
end

