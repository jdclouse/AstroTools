%% Monte Carlo Analysis
% Disturbance solar torque
analysis_set = 'MonteCarlo';
r = 35*pi/180;
offset_max = 1;
offset_min = -offset_max;
t = 0:1:Ts*3;

num_runs = 100;
y_int_mc_store = [];
y_lqr_mc_store = [];

for ii = 1:num_runs
    % Disturbed A
    fprintf('Run %d, ',ii)
    offset = offset_min + (offset_max-offset_min)*rand(1);
    fprintf('Offset %.4f\r',offset)
    dist_torque_max = offset*Fn_max;
    tmp = [zeros(1,4);
           dist_torque_max*B(2), zeros(1,3);
           zeros(1,4);
           dist_torque_max*B(4), zeros(1,3)];
    A_dist = A + tmp;
    
    A_OL_Aug_MC = [A_dist,zeros(4,1);-C, zeros(1)];
    A_CL_Aug_MC = [A_dist-B*K, -B*KI; -C, zeros(1)];
    A_Obs_Aug_MC = [A_OL_Aug_MC-B_OL_Aug*K_Aug,B_OL_Aug*K_Aug(1:4);
        zeros(4,5),A_dist-L*C];
    Int_sys_MC = ss(A_Obs_Aug_MC, B_Obs_Aug, C_Obs_AugFake, 0);
    
    A_Obs_LQR_MC = ...
        [A_OL_Aug_MC-B_OL_Aug*K_LQR,B_OL_Aug*K_LQR(1:4);...
         zeros(4,5),A_dist-L*C];
    LQR_system_MC = ss(A_Obs_LQR_MC, B_Obs_LQR, C_Obs_LQRFake, 0);
    
    y_int_mc = lsim(Int_sys_MC,repmat(r,1,length(t)),t);
    y_lqr_mc = lsim(LQR_system_MC,repmat(r,1,length(t)),t);
    
    y_int_mc_store(:,:,ii) = y_int_mc;
    y_lqr_mc_store(:,:,ii) = y_lqr_mc;
end
fprintf('\n')
analysis_set = 'LQRMonteCarlo';
plotSailSysResp( analysis_set,y_lqr_mc_store(:,1:5,:),t,K_LQR,r,Ts,3600*2 )
analysis_set = 'Ctrl1MonteCarlo';
plotSailSysResp( analysis_set,y_int_mc_store(:,1:5,:),t,K_Aug,r,Ts,3600*2 )
t = 0:0.01:Ts*3;