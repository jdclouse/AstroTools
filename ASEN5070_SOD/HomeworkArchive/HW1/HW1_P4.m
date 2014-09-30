%% HW1 Problem 4: Orbit Numerical Integration
fprintf('\n');
clearvars -except function_list pub_opt
close all

ode_opts = odeset('RelTol', 1e-12, 'AbsTol', 1e-20);
r = [-2436.45; -2436.45; 6891.037]; % km
v = [5.088611; -5.088611; 0.0]; % km/s
state = [r;v];

%Find the period
OE = cart2oe(state);
a = OE(1);
period = 2*pi*sqrt(a*a*a/3.986e5);

times = 0:20:period*2;

[T,X] = ode45(@two_body_state_dot, times, state, ode_opts);

%Get the magnitues for plotting
r_mag = zeros(1,length(times));
v_mag = zeros(1,length(times));
a_mag = zeros(1,length(times));

for i = 1:length(times)
    r_mag(i) = norm(X(i,1:3));
    v_mag(i) = norm(X(i,4:6));
    s_dot = two_body_state_dot(0,X(i,:));
    a_mag(i) = norm(s_dot(4:6));
end

%Plot the result
figHandle = figure;
set(figHandle, 'Position', [100, 100, 600, 800])
subplot(3,1,1);
plot(times, r_mag)
title('Position magnitude over two orbits')
ylabel('r (km)');
xlabel('time (s)');
subplot(3,1,2)
plot(times, v_mag)
title('Velocity magnitude over two orbits')
ylabel('v (km/s)');
xlabel('time (s)');
subplot(3,1,3)
plot(times, a_mag)
title('Acceleration magnitude over two orbits')
ylabel('a (km/s^2)');
xlabel('time (s)');

%% HW1 Problem 5: Orbit Numerical Integration Energy
% Why is the change in total specific energy not constant?
%
% The computational precision does not allow for exact calculations. The
% small error that results is built upon as the simulation runs.
KE = v_mag.*v_mag/2;
PE = -3.986e5./r_mag;

deltaE = KE + PE - (KE(1)+PE(1));
figure
plot(times, deltaE)

title('Change in energy over two numerically-integrated orbits')
ylabel('\DeltaEnergy');
xlabel('time (s)');
