%% Initialize
fprintf('\n');
clearvars -except function_list pub_opt
close all
% Bring in answers to compare
hw5_p1_answers

%% "Truth" solution
% Matlab's ode45 integrator was used, RelTol = 1e-12 and AbsTol = 1e-20

Xt0 = [1;0;0;1];
num_time_units = 10;
dt = 0.01; %TU
times = 0:dt:num_time_units; 

ode_opts = odeset('RelTol', 1e-13, 'AbsTol', 1e-20);
[T,Xout] = ode45(@hw5_deriv, times, Xt0, ode_opts);

% Record only at integer time units from t0
Xti = zeros(num_time_units+1, length(Xt0));
Xti(1,:) = Xt0';

for ii = 1:num_time_units
    Xti(ii+1,:) = Xout(ii/dt+1,:); %+1?
end
Xti(11,:)
Xti(11,:)'-X_10

%% Reference Trajectory
Xreft0 = Xt0 - [1e-6; -1e-6; 1e-6;1e-6];
STM = reshape(eye(4),16,1);

[T,Xout] = ode45(@hw5_deriv, times, [Xreft0; STM], ode_opts);

% Record only at integer time units from t0
Xrefti = zeros(num_time_units+1, length(Xt0));
STM_i=zeros(4,4, num_time_units+1); % 4x4xt
Xrefti(1,:) = Xreft0';
STM_i(:,:,1) = eye(4);

for ii = 1:num_time_units
    Xrefti(ii+1,1:4) = Xout(ii/dt+1,1:4); %+1?
    STM_i(:,:,ii+1) = reshape(Xout(ii/dt+1,5:20),4,4);
end
Xrefti(11,:)'-Xref_10
STM_i(:,:,end)-STM_10
Xti(11,:)'-Xrefti(11,:)'
STM_i(:,:,end)*(Xti(11,:)'-Xrefti(11,:)')%-STM_dX_10

%%
dim=2
J = [zeros(dim) eye(dim); -eye(dim) zeros(dim)];

inv(STM_i(:,:,end)) - -(J*STM_i(:,:,end)*J)'
inv(STM_i(:,:,end))*STM_i(:,:,end)