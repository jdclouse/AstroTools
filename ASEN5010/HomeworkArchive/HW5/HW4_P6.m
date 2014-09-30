%% Problem 6: S&J, Problem 8.8
fprintf('Problem 6: S&J, Problem 8.8\n');
clearvars -except function_list pub_opt
close all

%% xr(t)=0
state = [2;0];
K=1;
P=1;

delta_t = 0.01;
t_end = 60 - delta_t; % seconds

% Arrays for recording and plotting
t_mat = 0:delta_t:t_end+delta_t;
[rows, cols] = size(t_mat);
pos_mat = zeros(1,cols);
xr_mat = zeros(1,cols);
pos_mat(:,1) = state(1);
xr_mat(:,1) = 0;
idx = 2;

% anonymous fcn for CLD
cld = @(refxddot, k_gain, pos_err, p_gain, vel_err) ...
    (refxddot - k_gain*pos_err - p_gain*vel_err);
% euler
for t = 0:delta_t:t_end
    xr=0;
    xrdot=0;
    xrddot=0;
    dx=state(1)-xr;
    dxdot=state(2)-xrdot;
    
    state_deriv=[state(2);cld(xrddot,K,dx,P,dxdot)];

  state = state + delta_t*state_deriv;
  
  % Updating array
  pos_mat(:,idx) = state(1);
  xr_mat(:,idx) = xr;
  idx = idx + 1;
end

font_size=8;
figure
subplot(2,1,1)
plot(t_mat, [pos_mat;xr_mat]);
hold on
% plot(t_mat, xr_mat);
% title('MRP Propagation','FontSize',font_size)
xlabel('time(s)','FontSize',font_size)
ylabel('Position','FontSize',font_size)
legend('Position', 'Reference Position')
grid on
set(gca,'FontSize',font_size)
title('Controller performance, xr(t)=0')
fprintf('\n\n\n');
subplot(2,1,2)
plot(t_mat, [pos_mat-xr_mat]);
xlabel('time(s)','FontSize',font_size)
ylabel('Error','FontSize',font_size)
grid on
set(gca,'FontSize',font_size)
title('Position error')

%% xr(t)=0
state = [2;0];
K=1;
P=1;

delta_t = 0.01;
t_end = 60 - delta_t; % seconds

% Arrays for recording and plotting
t_mat = 0:delta_t:t_end+delta_t;
[rows, cols] = size(t_mat);
pos_mat = zeros(1,cols);
xr_mat = zeros(1,cols);
pos_mat(:,1) = state(1);
xr_mat(:,1) = 0;
idx = 2;

% euler
for t = 0:delta_t:t_end
    xr=sin(t/2);
    xrdot=cos(t/2)/2;
    xrddot=-sin(t/2)/4;
    dx=state(1)-xr;
    dxdot=state(2)-xrdot;
    
    state_deriv=[state(2);cld(xrddot,K,dx,P,dxdot)];

  state = state + delta_t*state_deriv;
  
  % Updating array
  pos_mat(:,idx) = state(1);
  xr_mat(:,idx) = xr;
  idx = idx + 1;
end

font_size=8;
figure
subplot(2,1,1)
plot(t_mat, pos_mat, t_mat, xr_mat);
% title('MRP Propagation','FontSize',font_size)
xlabel('time(s)','FontSize',font_size)
ylabel('Position','FontSize',font_size)
legend('Position', 'Reference Position')
grid on
set(gca,'FontSize',font_size)
title('Controller performance, xr(t)=sin(t/2)')
fprintf('\n\n\n');
subplot(2,1,2)
plot(t_mat, [pos_mat-xr_mat]);
xlabel('time(s)','FontSize',font_size)
ylabel('Error','FontSize',font_size)
grid on
set(gca,'FontSize',font_size)
title('Position error')
