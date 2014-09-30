%% Simulation with PD control
fprintf('\n');
clearvars -except function_list pub_opt
close all
MRP0=[0.8;0.1;-0.1];%rad

fprintf('MRP vector:\n');
MRP = MRP0; 
printVector(MRP, '');
omega_body0 = [0; 0; 0]; %rad/s
fprintf('Initial Body Rates:\n');
printVector(omega_body0, 'rad/s');
I_body=[10 0 0;0 20 0;0 0 30]; %kg*m2

%unit gains
P1=2*10/100;
K = P1*P1/10%0.004;
P = [P1 0 0; 0 sqrt(K*20) 0; 0 0 sqrt(K*30)]
cm_torque=[0;0;0];

delta_t = 0.01;
t_end = 600 - delta_t; % seconds

MRP = MRP0; 
omega_body = omega_body0; 
% Arrays for recording and plotting
t_mat = 0:delta_t:t_end+delta_t;
[rows, cols] = size(t_mat);
MRP_mat = zeros(3,cols);
omega_mat = zeros(3,cols);
EA_mat = zeros(3,cols);
mode_mat = zeros(3,cols);
CT_mat = zeros(3,cols);
MRP_mat(:,1) = MRP;
omega_mat(:,1) = omega_body;
mode_mat(:,1) = 0;
CT_MAT(:,1) = cm_torque;
idx = 2;

% RK4 integration
state = [MRP; omega_body];
control_int = 0;
linearization_mat = ...
    [zeros(3,3) 0.25*I_body; -K*inv(I_body) -inv(I_body)*P];
for t = 0:delta_t:t_end


  k1 = linearization_mat*state;
  k2 = linearization_mat*(state+delta_t*k1/2);
  k3 = linearization_mat*(state+delta_t*k2/2);
  k4 = linearization_mat*(state+delta_t*k3);

  state = state + delta_t/6*(k1 + 2*k2 + 2*k3 + k4);
      
  % Enforce |MRP| <= 1, switch to shadow set if needed
  if norm(state(1:3)) > 1
      state(1:3) = -state(1:3)/dot(state(1:3), state(1:3));
  end
  
  % Updating array
  
  MRP_mat(:,idx) = state(1:3);    
  omega_mat(:,idx) = state(4:6);
  idx = idx + 1;
end

font_size=8;
figure
plot(t_mat, MRP_mat);
mytitle = strcat('MRP Propagation - Linearized CLD');
title(mytitle,'FontSize',font_size)
xlabel('time(s)','FontSize',font_size)
ylabel('Element Magnitude','FontSize',font_size)
legend('\sigma_{1}', '\sigma_{2}', '\sigma_{3}')
grid on
set(gca,'FontSize',font_size)
fprintf('\n\n\n');

figure
plot(t_mat, omega_mat*180/pi);
mytitle = strcat('Body Rates - Linearized CLD');
title(mytitle,'FontSize',font_size)
xlabel('time(s)','FontSize',font_size)
ylabel('Element Magnitude (deg/s)','FontSize',font_size)
legend('\omega_{1}', '\omega_{2}', '\omega_{3}')
grid on
set(gca,'FontSize',font_size)
fprintf('\n\n\n');


%% Nonlinear CLD
% Arrays for recording and plotting
t_mat = 0:delta_t:t_end+delta_t;
[rows, cols] = size(t_mat);
MRP_mat = zeros(3,cols);
omega_mat = zeros(3,cols);
EA_mat = zeros(3,cols);
mode_mat = zeros(3,cols);
CT_mat = zeros(3,cols);
MRP_mat(:,1) = MRP;
omega_mat(:,1) = omega_body;
mode_mat(:,1) = 0;
CT_MAT(:,1) = cm_torque;
idx = 2;
% RK4 integration
state = [MRP; omega_body];
for t = 0:delta_t:t_end

  k1 = [derivMRP(state(1:3), state(4:6)); ...
      inv(I_body)*(-P*state(4:6)-K*state(1:3))];
  k2 = [derivMRP(state(1:3) + delta_t*k1(1:3)/2, ...
      state(4:6) + delta_t*k1(4:6)/2); ...
      inv(I_body)*(-P*(state(4:6)+delta_t*k1(4:6)/2)-K*(state(1:3)+delta_t*k1(1:3)/2))];
  k3 = [derivMRP(state(1:3) + delta_t*k2(1:3)/2, ...
      state(4:6) + delta_t*k2(4:6)/2); ...
      inv(I_body)*(-P*(state(4:6)+delta_t*k1(4:6)/2)-K*(state(1:3)+delta_t*k2(1:3)/2))];
  k4 = [derivMRP(state(1:3) + delta_t*k3(1:3), ...
      state(4:6) + delta_t*k3(4:6)); ...
     inv(I_body)*(-P*(state(4:6)+delta_t*k1(4:6)/2)-K*(state(1:3)+delta_t*k3(1:3)))];

  state = state + delta_t/6*(k1 + 2*k2 + 2*k3 + k4);
      
  % Enforce |MRP| <= 1, switch to shadow set if needed
  if norm(state(1:3)) > 1
      state(1:3) = -state(1:3)/dot(state(1:3), state(1:3));
  end
  
  % Updating array
  
  MRP_mat(:,idx) = state(1:3);    
  omega_mat(:,idx) = state(4:6);
%   CT_mat(:,idx)=control_torque;
  idx = idx + 1;
end

font_size=8;
figure
plot(t_mat, MRP_mat);
mytitle = strcat('MRP Propagation - Nonlinear CLD');
title(mytitle,'FontSize',font_size)
xlabel('time(s)','FontSize',font_size)
ylabel('Element Magnitude','FontSize',font_size)
legend('\sigma_{1}', '\sigma_{2}', '\sigma_{3}')
grid on
set(gca,'FontSize',font_size)
fprintf('\n\n\n');

figure
plot(t_mat, omega_mat*180/pi);
mytitle = strcat('Body Rates - Nonlinear CLD');
title(mytitle,'FontSize',font_size)
xlabel('time(s)','FontSize',font_size)
ylabel('Element Magnitude (deg/s)','FontSize',font_size)
legend('\omega_{1}', '\omega_{2}', '\omega_{3}')
grid on
set(gca,'FontSize',font_size)
fprintf('\n\n\n');
