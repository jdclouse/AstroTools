%% Monte Carlo Analysis
%% Disturbance solar torque
analysis_set = 'MonteCarlo';
r = 35*pi/180;
offset_max = 0.1;
offset_min = -offset_max;
t = 0:1:Ts*3;

num_runs = 100;
y_int_mc_store = [];
y_lqr_mc_store = [];

for ii = 1:num_runs
    fprintf('Run %d, ',ii)
    offset = offset_min + (offset_max-offset_min)*rand(1);
    fprintf('Offset %.4f\r',offset)
    dist_torque_max = offset*Fn_max;
    A_dist = A;
    % The disturbance will manifest itself as an input torque
    B_dist = [0; 1/a_dd_LHS_1; 0; -BIG/d_dd_LHS;zeros(5,1)];
        
    A_OL_Aug_MC = [A_dist,zeros(4,1);-C, zeros(1)];
    A_CL_Aug_MC = [A_dist-B*K, -B*KI; -C, zeros(1)];
    A_Obs_Aug_MC = [A_OL_Aug_MC-B_OL_Aug*K_Aug,B_OL_Aug*K_Aug(1:4);
        zeros(4,5),A_dist-L*C];
    Int_sys_MC = ss(A_Obs_Aug_MC, [B_Obs_Aug,B_dist], C_Obs_AugFake, 0);
    
    A_Obs_LQR_MC = ...
        [A_OL_Aug_MC-B_OL_Aug*K_LQR,B_OL_Aug*K_LQR(1:4);...
         zeros(4,5),A_dist-L*C];
    LQR_system_MC = ss(A_Obs_LQR_MC, [B_Obs_LQR,B_dist], C_Obs_LQRFake, 0);
    
    y_int_mc = lsim(Int_sys_MC,repmat([r;dist_torque_max],1,length(t)),t);
    y_lqr_mc = lsim(LQR_system_MC,repmat([r;dist_torque_max],1,length(t)),t);
%     y_int_mc = lsim(Int_sys_MC,repmat(r,1,length(t)),t);
%     y_lqr_mc = lsim(LQR_system_MC,repmat(r,1,length(t)),t);
    
    y_int_mc_store(:,:,ii) = y_int_mc;
    y_lqr_mc_store(:,:,ii) = y_lqr_mc;
end
fprintf('\n')
analysis_set = 'LQRMonteCarlo';
plotSailSysResp( analysis_set,y_lqr_mc_store(:,1:5,:),t,K_LQR,r,Ts,3600*2,...
    title_plots)
analysis_set = 'Ctrl1MonteCarlo';
plotSailSysResp( analysis_set,y_int_mc_store(:,1:5,:),t,K_Aug,r,Ts,3600*2,...
    title_plots)
t = 0:0.01:Ts*3;

%% LQR update - Make it meet the gimbal bounds with the disturbance
Q_wts = [1,1,11000,1,1];
Q_wts = Q_wts/sum(Q_wts);
state_max = [pi/2, 0.01, pi/6, 0.01, 0.01];
Q = diag(Q_wts.*Q_wts./(state_max.*state_max));
rho_R = 1000;
u_max = 100;
R = rho_R/u_max;
[K_LQR, W, E] = lqr(A_OL_Aug,B_OL_Aug,Q,R);

A_Obs_LQR = [A_OL_Aug-B_OL_Aug*K_LQR,B_OL_Aug*K_LQR(1:4);zeros(4,5),A-L*C];
B_Obs_LQR = [zeros(size(B));1;zeros(length(L),1)];
B_dist = [0; 1/a_dd_LHS_1; 0; -BIG/d_dd_LHS;zeros(5,1)];
C_Obs_LQR = [C, 0, zeros(1,length(L))];
C_Obs_LQRFake = [eye(9)];
y_lqr_mc_store = [];
t = 0:1:Ts*3;
for ii = 1:num_runs
    fprintf('Run %d, ',ii)
    offset = offset_min + (offset_max-offset_min)*rand(1);
    fprintf('Offset %.4f\r',offset)
    dist_torque_max = offset*Fn_max;
    
    LQR_system = ss(A_Obs_LQR, [B_Obs_LQR,B_dist], C_Obs_LQRFake, 0);
    y_lqr_mc = lsim(LQR_system,repmat([r;dist_torque_max],1,length(t)),t);
    
    y_lqr_mc_store(:,:,ii) = y_lqr_mc;
end
analysis_set = 'LQR2MonteCarlo';
plotSailSysResp( analysis_set,y_lqr_mc_store(:,1:5,:),t,K_LQR,r,Ts,3600*2,...
    title_plots)
t = 0:0.01:Ts*3;