%% HW3 Problem 3
fprintf('\n');
clearvars -except function_list hw_pub toolsPath 
close all
CelestialConstants; % import useful constants

n = 15.5918272 / day2sec *2*pi;
M = 231.8021 * pi/180;
w = 106.9025 * pi/180;
e = 0.0008148;
RAAN = 342.1053 * pi/180;
i = 51.6396 * pi/180;

M2 = M + 3600*n;
if M2 > 2*pi
    M2 = M2 - 2*pi;
end

f = E2f(M2E(M2,e),e);
a = (Earth.mu/n^2)^(1/3);

[r, v ] = OE2cart( a,e,i,RAAN,w,f,Earth.mu);

fprintf('r,v for ISS one hour after epoch:\n');
fprintf('r_ECI = %.2f\n',r(1));
fprintf('        %.2f km\n',r(2));
fprintf('        %.2f\n\n',r(3));
fprintf('v_ECI = %.3f\n',v(1));
fprintf('        %.3f km/s\n',v(2));
fprintf('        %.3f\n',v(3));