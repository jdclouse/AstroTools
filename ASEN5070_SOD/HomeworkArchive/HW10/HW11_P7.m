%% Markov
%% Initialize
clearvars -except function_list pub_opt P_joseph_store
global function_list;
function_list = {};
close all

obs_data = load('hw11.dat');

T = 10;
truth = sin(obs_data(:,1) * 2*pi/T);
figure
hold on
plot(obs_data(:,1),obs_data(:,2))
plot(obs_data(:,1),truth,'r','LineWidth',3)

eta0_ap = 0;
P0_ap = 1;
R = 1;
Q = 1;
num_obs = length(obs_data(:,1));
eta_est_store = zeros(num_obs,1);

eta_est = eta0_ap;
P = P0_ap;

tc_vec = [.1 1 3 5];
sig_vec = [.5 4 2 1];

best_RMS = -1;
for idx = 1:length(tc_vec)
time_const = tc_vec(idx);
beta = 1/time_const;
sigma = sig_vec(idx);
eta_est_store_inner = zeros(num_obs,1);
RMS_accum = 0;
for ii = 1:num_obs
    % STM
    if ii == 1 %measurement at t = 0
        m = 1;
    else
        m = exp(-beta*(obs_data(ii,1)-obs_data(ii-1,1)));
    end
    STM = m;
    
    % Time Update
    eta_ap = STM*eta_est;
    gamma = sqrt(sigma*sigma/2/beta*(1-m*m));
    P_ap = STM*P*STM + gamma*Q*gamma;
    
    % Kalman gain
    K = P_ap/(P_ap+1); % valid for this 1D case
    
    % Measurement Update
    Y = obs_data(ii,2);
    eta_est = eta_ap +K*(Y-eta_ap); %H~ == 1 in this case.
    P = K;
    
    eta_est_store_inner(ii) = eta_est;
    RMS_accum = RMS_accum + (truth(ii)-eta_est)*(truth(ii)-eta_est);
    
end
RMS = sqrt(RMS_accum/num_obs);
fprintf(sprintf('RMS for tau=%f, sigma=%f: %f\n',time_const, sigma, RMS));
if best_RMS == -1
    best_RMS = RMS;
    eta_est_store = eta_est_store_inner;
else
    if RMS < best_RMS
        eta_est_store = eta_est_store_inner;
    end
    best_RMS = min(best_RMS,RMS);
end
end
plot(obs_data(:,1),eta_est_store,'g','LineWidth',2)
legend('Observations', 'Truth', sprintf('Gauss-Markov, RMS = %f',best_RMS))
% figure
% hist(obs_data(ii,2)-eta_est_store);