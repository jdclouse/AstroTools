%% HW5 Problem 4
fprintf('\n');
clearvars -except function_list hw_pub toolsPath 
close all
CelestialConstants; % import useful constants

t = 300; % sec
n = sqrt(Earth.mu/(400+Earth.R)^3);
rel_state = [200;300;50;0.1;0;-0.1]; %m
final_rel_state = [0;0;0;0;0;0]; %m
figure
color = 'b';

for t = [300, 30*60]
% Can take the state transition matrix, divide it into sub-matrices
% and solve for the required initial velocities.
% Subtract the contribution from initial pos to the final pos, from the 
% final pos. 
% Multiply by the inverse of the contribution of the initial vel to the
% final pos.
fprintf('t = %d sec:\n',t);
STM = CWHillSTM(n,t);
req_vel_init = inv(STM(1:3, 4:6))*...
    (final_rel_state(1:3)-STM(1:3, 1:3)*rel_state(1:3));
dv1 = norm(req_vel_init - rel_state(4:6))
final_state_preburn = STM*[rel_state(1:3);req_vel_init];
dv2 = norm(-final_state_preburn(4:6))
dv_tot = dv1 + dv2
for ii = 1:t
    plot_state = CWHillSTM(n,ii)*[rel_state(1:3);req_vel_init];
    x(ii) = plot_state(1);
    y(ii) = plot_state(2);
end
plot(y,x,color)
axis equal
hold on
color = 'r';
end
xlabel('In-Track')
ylabel('Radial')
legend('5 minutes', '30 minutes')