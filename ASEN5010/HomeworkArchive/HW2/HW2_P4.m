%% Problem 4: S&J, Problem 3.12
fprintf('Problem 4: S&J, Problem 3.12')
clearvars -except function_list pub_opt
close all
fprintf('Initial Euler Angles (3-2-1):\n')
euler_angles = [40; 30; 80]; % degrees
printVector(euler_angles, 'degrees');

delta_t = 0.01;
t_end = 60 - delta_t; % seconds

% Arrays for recording and plotting
t_mat = 0:delta_t:t_end+delta_t;
[rows, cols] = size(t_mat);
EA_mat = zeros(3,cols);
EA_mat(:,1) = euler_angles;
idx = 2;

% Euler integration
for t = 0:delta_t:t_end
  
  w_body_frame = [sin(0.1*t); 0.01; cos(0.1*t)] * 20 * pi/180; % rad/s
  euler_angles_dot = BmatEuler('321', euler_angles*pi/180)*w_body_frame; 
      %rad/s
  euler_angles = euler_angles + euler_angles_dot * delta_t * 180/pi; 
      % degrees, position at t_(n+1)

  % Updating array
  EA_mat(:,idx) = euler_angles;
  
  idx = idx + 1;
end

% The answer
fprintf('Euler angles after 60 seconds:\n')
printVector( euler_angles, 'degrees')

plot(t_mat, EA_mat);
title('3-2-1 Euler Angles')
xlabel('time(s)')
ylabel('degrees')
legend('yaw', 'pitch', 'roll')
grid on
fprintf('\n\n\n');