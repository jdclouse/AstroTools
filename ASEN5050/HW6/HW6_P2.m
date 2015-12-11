%% HW6 Problem 2
fprintf('\n');
% clearvars -except function_list hw_pub toolsPath 
% hold on
close all
% CelestialConstants; % import useful constants

lat = 40.01 * pi/180; 
lon = 254.83 * pi/180;
alt = 1.615;
topo_store = zeros(3,length(time_vec));

% Store off all the topographic coord vectors from the r_ecef info in the
% last problem
for ii = 1:length(r_ecef_store(1,:))
    topo_store(:,ii) = ecef2topo( r_ecef_store(:,ii), lat, lon, alt);
end

% Get the first pass indices. Elevation is the second element of each
% vector. Look for them being > 0.
first_idx = find(topo_store(2,:) > 0,1);
last_idx =  find(topo_store(2,first_idx:end) <= 0,1)+first_idx-2;

% Warning, this part assumes dt = 1 minute, and that the first pass happens
% in the first hour of the scenario.
first_time = (first_idx - 1); % minutes
last_time = (last_idx - 1); % minutes
fprintf(['Pass begins: 10-11-2014 0' num2str(1) ':' ...
    num2str(first_time) ':00 GMT\n']);
fprintf(['Pass ends:   10-11-2014 0' num2str(1) ':' ...
    num2str(last_time) ':00 GMT\n']);
fprintf(['Pass duration is approx ' num2str(last_time-first_time) ...
    ' minutes\n'])
fprintf(['max elevation is ' ...
    num2str(max(topo_store(2,first_idx:last_idx)*180/pi)) ...
    ' degrees\n'])

% A plot with North being up on the plot, clockwise angles
figure('OuterPosition', [0 50 hw_pub.figWidth hw_pub.figHeight])
plotVisZenith(topo_store(1,first_idx:last_idx), ...
    topo_store(2,first_idx:last_idx))
title('ISS Visibility from Boulder, CO (nice plot)')

% regular ol' polar plot
figure('OuterPosition', [0 50 hw_pub.figWidth hw_pub.figHeight])
polar(topo_store(1,first_idx:last_idx), ...
    90-topo_store(2,first_idx:last_idx)*180/pi)
xlabel('Azimuth (deg), Zenith (deg)')
title('ISS Visibility from Boulder, CO')