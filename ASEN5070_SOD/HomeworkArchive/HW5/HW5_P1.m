%% HW5 Problem 1
%% Initialize
fprintf('\n');
clearvars -except function_list pub_opt
close all
% Bring in answers to compare
hw5_p1_answers

%% a) "Truth" solution
% Matlab's ode45 integrator was used, RelTol = 1e-12 and AbsTol = 1e-20.
% Time interval of 0.01 Time Units used for integration.

Xt0 = [1;0;0;1];
num_time_units = 100;
dt = 0.01; %TU
times = 0:dt:num_time_units; 

ode_opts = odeset('RelTol', 1e-13, 'AbsTol', 1e-20);
[T,Xout] = ode45(@hw5_deriv, times, Xt0, ode_opts);

% Record only at integer time units from t0
num_record_step=10;
Xti = zeros(num_record_step+1, length(Xt0));
Xti(1,:) = Xt0';

for ii = 1:num_record_step
    Xti(ii+1,:) = Xout(num_record_step*ii/dt+1,:); %+1?
end

%Comparison
% fprintf('Nominal diffs:\n');
% Xti(2,:)'-X_10
% Xti(11,:)'-X_100

fprintf('Truth trajectory, t=10 TU:\n')
for ii = 1:4
    fprintf('%.9f\n', Xti(2,ii))
end
fprintf('\nTruth trajectory, t=100 TU:\n')
for ii = 1:4
    fprintf('%.9f\n', Xti(11,ii))
end
fprintf('\n')

%% b) Reference Trajectory
Xreft0 = Xt0 - [1e-6; -1e-6; 1e-6;1e-6];
STM = reshape(eye(4),16,1);

[T,Xout] = ode45(@hw5_deriv, times, [Xreft0; STM], ode_opts);

% Record only at integer time units from t0
Xrefti = zeros(num_record_step+1, length(Xt0));
STM_i=zeros(4,4, num_record_step+1); % 4x4xt
Xrefti(1,:) = Xreft0';
STM_i(:,:,1) = eye(4);

for ii = 1:num_record_step
    Xrefti(ii+1,1:4) = Xout(num_record_step*ii/dt+1,1:4); %+1?
    STM_i(:,:,ii+1) = reshape(Xout(num_record_step*ii/dt+1,5:20),4,4);
end

%Comparison
% fprintf('Ref trajectory, STM diffs:\n');
% Xrefti(2,:)'-Xref_10
% STM_i(:,:,2)-STM_10
% Xrefti(11,:)'-Xref_100
% STM_i(:,:,end)-STM_100
% Xti(2,:)'-Xrefti(2,:)'-dX_10
% Xti(11,:)'-Xrefti(11,:)'-dX_100
% STM_i(:,:,2)*(Xti(1,:)'-Xrefti(1,:)')-STM_dX_10
% STM_i(:,:,end)*(Xti(1,:)'-Xrefti(1,:)')-STM_dX_100
fprintf('Reference trajectory, t=10 TU:\n')
for ii = 1:4
    fprintf('%.9f\n', Xrefti(2,ii))
end
fprintf('\nReference trajectory, t=100 TU:\n')
for ii = 1:4
    fprintf('%.9f\n', Xrefti(11,ii))
end
fprintf('\n')
%% c) Show STM is symplectic
dim=2;
J = [zeros(dim) eye(dim); -eye(dim) zeros(dim)];

inv_STM = -(J*STM_i(:,:,end)*J)';
fprintf('STM inverse, t=100 TU:\n')
disp(inv_STM)

prod = STM_i(:,:,end)*inv_STM;
fprintf('\nSTM*inv(STM), t=100 TU:\n')
for ii = 1:4
fprintf('%.9f %.9f %.9f %.9f\n', prod(ii,1), prod(ii,2), prod(ii,3), ...
    prod(ii,4));
end

%% d) Calculate perturbation vector
% The different methods of calculating dX are pretty small (<0.1%), different due
% to the numerical propagation of the truth, reference, and STM.
dX_method1 = Xti(11,:)'-Xrefti(11,:)'
dX_method2 = STM_i(:,:,end)*(Xti(1,:)'-Xrefti(1,:)')
dX_diff = dX_method1 - dX_method2
