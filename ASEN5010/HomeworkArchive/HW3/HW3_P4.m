%% Problem 4: S&J, Problem 3.28
fprintf('Problem 4: S&J, Problem 3.28\n');
clearvars -except function_list pub_opt
close all
fprintf('MRP vector:\n');
MRP = [0; 0; 0]; 
printVector(MRP, '');
omega_body = [1; 0.5; -0.7];

delta_t = 0.01;
t_end = 5 - delta_t; % seconds

% Arrays for recording and plotting
t_mat = 0:delta_t:t_end+delta_t;
[rows, cols] = size(t_mat);
MRP_mat = zeros(3,cols);
MRP_mat(:,1) = MRP;
idx = 2;

% Checking work... I know the 321 EA works.
phi = 2*acos((1-dot(MRP,MRP))/(1+dot(MRP,MRP)));
PRV = phi * MRP/tan(phi/4);
% euler_angles = DCM2Euler('321', PRV2DCM(PRV))
euler_angles = MRP
EA_mat = zeros(3,cols);
EA_mat(:,1) = euler_angles;

% Euler integration
for t = 0:delta_t:t_end
  MRP_dot = derivMRP(MRP, omega_body); 
      %rad/s
  MRP = MRP + MRP_dot * delta_t; 
      % degrees, position at t_(n+1)
      
  % Enforce |MRP| <= 1, switch to shadow set if needed
  if norm(MRP) > 1
      MRP = -MRP/dot(MRP, MRP);
  end
  
  euler_angles_dot = BmatEuler('321', euler_angles)*omega_body; 
      %rad/s
  euler_angles = euler_angles + euler_angles_dot * delta_t;
  EA_mat(:,1) = euler_angles;
  
  % Updating array
  MRP_mat(:,idx) = MRP;
  
  idx = idx + 1;
end
% The answer
fprintf('MRP elements after 5 seconds:\n')
printVector( MRP, '')
euler_angles
phi = 2*acos((1-dot(MRP,MRP))/(1+dot(MRP,MRP)));
PRV = phi * MRP/tan(phi/4);
euler_angles_from_MRP = DCM2Euler('321', PRV2DCM(PRV))
euler_angles-euler_angles_from_MRP

plot(t_mat, MRP_mat);
title('MRP Elements')
xlabel('time(s)')
ylabel('Element Magnitude')
legend('\sigma_{1}', '\sigma_{2}', '\sigma_{3}')
grid on
fprintf('\n\n\n');