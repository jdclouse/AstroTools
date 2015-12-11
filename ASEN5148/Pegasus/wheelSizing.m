%% Wheel Sizing
I = 8e-4; %Surrey SP-10
w_min = 500*2*pi/60; %r/s
w_max = 5000*2*pi/60; %r/s

H_min = I*w_min
H_max = I*w_max

% max disturbance over orbit
max_mom_change = 1e-5*2*pi*sqrt((6978^3)/3.986e5)*0.707
% Time to completely damp
mag_moment = 0.001;
damp_time = H_max/mag_moment