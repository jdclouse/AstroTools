%% HW1 Problem 4: Ground station observations
fprintf('\n');
clearvars -except function_list hw_pub toolsPath 
close all

v = 50; %m/s
h = 100; %m
x0 = 250; %m

x = linspace(-x0, x0, x0*10);

range = sqrt(h*h + x.*x);
range_rate = x./range*v;
zenith = atan2(x,h);

h2 = 50; %m
range2 = sqrt(h2*h2 + x.*x);
range_rate2 = x./range2*v;
zenith2 = atan2(x,h2);

figure('OuterPosition', [0 50 hw_pub.figWidth hw_pub.figHeight])
subplot(3,1,1)
plot(x, range)
hold on
plot(x, range2, 'g')
title('Range')
xlabel('x position (m)')
ylabel('range (m)')
legend('h = 100 m', 'h = 50 m')
subplot(3,1,2)
plot(x, range_rate)
hold on
plot(x, range_rate2, 'g')
title('Range Rate')
xlabel('x position (m)')
ylabel('range rate (m/s)')
legend('h = 100 m', 'h = 50 m')
subplot(3,1,3)
plot(x, zenith*180/pi)
hold on
plot(x, zenith2*180/pi, 'g')
title('Zenith Angle')
xlabel('x position (m)')
ylabel('zenith (deg)')
legend('h = 100 m', 'h = 50 m')

fprintf(['range:\n',...
    '\tthe slope is steep away from the ground station, making the\n',...
    '\tmeasurement more sensitive to changes in x (more accurate)\n.',...
    '\tHowever, it loses accuracy directly above. There are also not\n',...
    '\tunique x/range pairings, so another method must be used to \n',...
    '\tdetermine which side of the transmitter the receiver is on.\n',...
    '\tLower values of h make the range measurement more accurate \n',...
    '\tcloser to the transmitter\n',...
    'range rate:\n',...
    '\tThe slope is steep only near the transmitter, so it is not\n',...
    '\taccurate at long distances.\n',...
    '\tIt has unique x/range rate pairings, so you can tell which side\n',...
    '\tof the transmitter the receiver is.\n',...
    '\tLower values of h make the range rate measurement less accurate\n',...
    '\tuntil the reciever is closest to the transmitter.\n',...
    'Zenith angle:\n',...
    '\tThe slope is steep only near the transmitter, so it is not\n',...
    '\taccurate at long distances.\n',...
    '\tIt has unique x/range rate pairings, so you can tell which side\n',...
    '\tof the transmitter the receiver is.\n',...
    '\tLower values of h make the range rate measurement less accurate\n',...
    '\tuntil the reciever is closest to the transmitter.\n'])

%% Hyperbolic multilateration
figure
hold on
color_opts = ['b';'g'; 'r'];
counter = 1;
for ground_station_spacing = [50 100]; 
    x1 = -ground_station_spacing/2; %m
    x2 = ground_station_spacing/2; %m
    r1 = sqrt((x1-x).*(x1-x) + h*h);
    r2 = sqrt((x2-x).*(x2-x) + h*h);
    d21 = r2 - r1;
    d21_dot = -(x2-x)./r2*v - -(x1-x)./r1*v;
    subplot(2,1,1)
    hold on
    plot(x, d21, color_opts(counter))
    title('Range Difference')
    xlabel('x position (m)')
    ylabel('\Deltarange (m)')
    legend('50 m tx distance', '100 m tx distance')
    subplot(2,1,2)
    hold on
    plot(x, d21_dot, color_opts(counter))
    title('Range Difference Rate')
    xlabel('x position (m)')
    ylabel('\Deltarange rate (m/s)')
    legend('50 m tx distance', '100 m tx distance')
    counter = counter + 1;
end
fprintf(['delta range:\n',...
    '\tthe slope is steep closer to the ground station, making the\n',...
    '\tmeasurement more sensitive to changes in x when it approaches \n',...
    '\tthe first transmitter.\n',...
    '\tThere are also unique x/delta-range pairings, so it is known\n',...
    '\twhere the receiver is wrt the transmitters\n',...
    '\tLarger transmitter distances result in more accuracy both between\n',...
    '\tand slightly beyond the transmitters\n',...
    'delta range rate:\n',...
    '\tThe slope is steep just as it approaches the transmitters.\n',...
    '\tIt does not have unique x/delta range rate pairings, so you cannot\n',...
    '\ttell which side of the midpoint the receiver is located.\n',...
    'Delta Range measurement seems most sufficient in this case\n'])