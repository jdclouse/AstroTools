clearvars -except function_list pub_opt
close all

stat_od_proj_init
ObsData = load('ObsData.txt');
A_given = load('BatchA.mat');

consts.Re = Re;
consts.area = drag.A;
consts.rho = compute_density(ri);
consts.theta_dot = theta_dot;
consts.m = drag.m;
consts.state_len = 18;

A = stat_od_proj_A(state, consts);

relDiffA = abs((A(1:6,1:9) -A_given.A)./A_given.A);


H_given = load('BatchHtilda.mat');

consts.t = 0;
consts.site = 2;
H_tilda = stat_od_proj_H_tilda(state, consts);
relDiffHtilda = abs((H_tilda -H_given.Htildes(:,:,1))./H_given.Htildes(:,:,1));

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
relDiffSTM20(1:6, 1:9)

STM_end = eye(consts.state_len);
STM_end(1:important_block(1),1:important_block(2)) = ...
    reshape(X(end,consts.state_len+1:end), ...
    important_block(1), important_block(2));
relDiffSTM20 = abs((STM_end-Phi_given.Phi_t18340)./Phi_given.Phi_t18340);
relDiffSTM20(1:6, 1:9)

%%
num_obs = length(ObsData(:,1));
RMS_r = 0;
RMS_rr = 0;
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
    
    RMS_r = RMS_r + (ObsData(ii,3)-r_comp)*(ObsData(ii,3)-r_comp);
    RMS_rr = RMS_rr + (ObsData(ii,4)-rr_comp)*(ObsData(ii,4)-rr_comp);
end
RMS_r = sqrt(RMS_r/num_obs)
RMS_rr = sqrt(RMS_r/num_obs)