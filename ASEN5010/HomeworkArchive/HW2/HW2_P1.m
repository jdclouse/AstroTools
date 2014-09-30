%% Problem 1: S&J, Problem 3.11
fprintf('Problem 1: S&J, Problem 3.11\n');
clearvars -except function_list pub_opt
close all
fprintf('Initial Euler Angles (3-2-1):\n');
euler_angles = [10; -15; 20]; % degrees
printVector(euler_angles, 'degrees');
fprintf('Initial Euler Angle Rates (3-2-1):\n');
euler_angle_rates = [2; 1; 0;]; % deg/s
printVector(euler_angle_rates, 'degrees/s');

%N->B rotation [BN]
BN = Euler2DCM('321', euler_angles * pi/180); 

%Find invers B matrix
w_body_B_frame = BinvEuler('321', euler_angles*pi/180) * euler_angle_rates*pi/180; % rad/s
fprintf('Omega body frame:\n');
printVector(w_body_B_frame * 180/pi, 'deg/s'); % deg/s

% w_body_N_frame = [NB]*w_body_B_frame
w_body_N_frame = inv(BN)*w_body_B_frame; 
fprintf('Omega N-frame:\n');
printVector(w_body_N_frame * 180/pi, 'deg/s'); % deg/s
fprintf('\n\n\n');