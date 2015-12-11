%% HW2 Problem 4
fprintf('\n');
clearvars -except function_list hw_pub toolsPath 
close all
CelestialConstants; % import useful constants

r = [-5650;-2650;2850]; %km
v = [2.415;-7.032;-1.796]; %km/s

[a,e,i,RAAN,w,f] = cart2OE(r,v,Earth.mu);

rp = a*(1-e);
hp = rp - Earth.R;
ra = a*(1+e);
ha = ra - Earth.R;

specific_energy = norm(v)*norm(v)/2-Earth.mu/norm(r);

h = cross(r,v);

fpa = atan2(e*sin(f),1+e*cos(f));

fprintf('a. Semi-major axis: %.2f km\n', a)
fprintf('b. Eccentricity: %.2f\n', e)
fprintf('c. Perigee altitude: %.2f km\n', hp)
fprintf('   Apogee altitude: %.2f km\n', ha)
fprintf('d. Specific energy: %.2f J/kg\n', specific_energy)
fprintf('e. Angular momentum: %.2f km2/s\n', h(1))
fprintf('                     %.2f km2/s\n', h(2))
fprintf('                     %.2f km2/s\n', h(3))
fprintf('f. Inclination: %.2f deg\n', i*180/pi)
fprintf('b. Flight path angle: %.2f deg\n', fpa*180/pi)