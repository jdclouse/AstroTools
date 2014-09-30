%% Initialize
fprintf('\n');
clearvars -except function_list pub_opt
close all

%% "Truth" solution
% Matlab's ode45 integrator was used, RelTol = 1e-12 and AbsTol = 1e-20

Xt0 = [1;0;0;1];
num_time_units = 10;
dt = 0.1; %TU
times = 0:dt:num_time_units; 

ode_opts = odeset('RelTol', 1e-12, 'AbsTol', 1e-20);
[T,Xout] = ode45(@hw5_deriv, times, Xt0, ode_opts);

% Record only at integer time units from t0
Xti = zeros(num_time_units+1, length(Xt0));
Xti(1,:) = Xt0';

for ii = 1:num_time_units
    Xti(ii+1,:) = Xout(ii/dt+1,:); %+1?
end

%% Reference Trajectory
Xreft0 = Xt0 - [1e-6; -1e-6; 1e-6;1e-6];
