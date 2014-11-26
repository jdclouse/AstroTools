%% HW 6: Setting Up the Term Project
%% Initialize
clearvars -except function_list pub_opt
ComputeA;
ComputeH_tilda;
clearvars -except function_list pub_opt
close all

stat_od_proj_init
ObsData = load('ObsData.txt');

%% Problem 1: Derive A, Htilda matrices
consts.Re = Re;
consts.area = drag.A;
consts.rho = compute_density(ri);
consts.theta_dot = theta_dot;
consts.m = drag.m;
consts.state_len = 18;

A_given = load('BatchA.mat');
A = stat_od_proj_A(state, consts);

relDiffA = abs((A(1:6,1:9) -A_given.A)./A_given.A);
fprintf('For A, elements that were exactly the same as the solution \n')
fprintf('cause Infs to show up in the log10 matrix. I replaced them\n')
fprintf('with zeros so that they were distinct from the other elements\n')
fprintf('in the histogram.\n')
hist_mat = log10(relDiffA);
hist_mat(hist_mat == -Inf) = 0;
hist(reshape(hist_mat,6*9,1));
title('A Matrix relative diff histogram.')
xlabel('Exponent')
ylabel('Num Elements')

H_given = load('BatchHtilda.mat');

consts.t = 0;
consts.site = 2;
H_tilda = stat_od_proj_H_tilda(state, consts);
relDiffHtilda = ...
    abs((H_tilda -H_given.Htildes(:,:,1))./H_given.Htildes(:,:,1));
relDiffHtilda

%% Problem 2: Integrate position, velocity, STM
dt = 0.1;
times = 0:dt:18340;
ode_opts = odeset('RelTol', 1e-12, 'AbsTol', 1e-20);
[T,X] = ode45(@two_body_state_dot, times, state, ode_opts, propagator_opts);

% Test the STMs against given ones
Phi_given = load('BatchPhi.mat');
STM_20 = eye(consts.state_len);
STM_20(1:important_block(1),1:important_block(2)) = ...
    reshape(X(20/dt+1,consts.state_len+1:end), ...
    important_block(1), important_block(2));
relDiffSTM20 = abs((STM_20-Phi_given.Phi_t20)./Phi_given.Phi_t20);
hist_mat = log10(relDiffSTM20);
hist_mat(hist_mat == -Inf) = 0;
figure
hist(reshape((hist_mat(1:6, 1:9)),6*9,1))
title('I(20,0) Matrix relative diff histogram.')
xlabel('Exponent')
ylabel('Num Elements')

STM_end = eye(consts.state_len);
STM_end(1:important_block(1),1:important_block(2)) = ...
    reshape(X(end,consts.state_len+1:end), ...
    important_block(1), important_block(2));
relDiffSTMend = abs((STM_end-Phi_given.Phi_t18340)./Phi_given.Phi_t18340);
% relDiffSTM20(1:6, 1:9)
figure
hist(reshape(log10(relDiffSTMend(1:6, 1:9)),6*9,1))
title('I(18340,0) Matrix relative diff histogram.')
xlabel('Exponent')
ylabel('Num Elements')

%% 
num_obs = length(ObsData(:,1));
r_residual_mat = [];
rr_residual_mat = [];
for ii = 1:num_obs
    site_num = 0;
    for jj = 1:3
        if ObsData(ii, 2) == site(jj).id
            site_num = jj;
            break
        end
    end
    t_obs = ObsData(ii,1);
    ostate = X(T(:,1)==t_obs,1:6);
    
    r_comp = compute_range_ECFsite(ostate(1:3),...
        site(site_num).r,theta_dot*t_obs);
    rr_comp = compute_range_rate_ECFsite(ostate(1:6),...
        site(site_num).r,theta_dot*t_obs, theta_dot);
    
    r_residual_mat(ii) = (ObsData(ii,3)-r_comp);
    rr_residual_mat(ii) = (ObsData(ii,4)-rr_comp);

end
figure
subplot(2,1,1)
plot(ObsData(:,1),r_residual_mat,'.')
ylabel('Range Residual')
title('Measurement Residuals')
subplot(2,1,2)
plot(ObsData(:,1),rr_residual_mat,'.')
ylabel('Range Rate Residual')
xlabel('Time Since Epoch (s)')
RMS_r_residual = sqrt(sum(r_residual_mat.*r_residual_mat)/num_obs);
RMS_rr_residual = sqrt(sum(rr_residual_mat.*rr_residual_mat)/num_obs);
fprintf('RMS range residual: %f m \n', RMS_r_residual)
fprintf('RMS range rate residual: %f m/s \n', RMS_rr_residual)