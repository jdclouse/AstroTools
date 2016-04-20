
num_obs = length(ObsMassaged(:,1));
% obs_diffs = ObsMassaged(2:end,2) - ObsMassaged(1:end-1,2);
% sum(obs_diffs>60)
% % The values in here are the indices where an observation segment ends.
% obs_seg_idx = [find(obs_diffs>60); num_obs];
obs_diffs = ObsMassaged(2:end,1) - ObsMassaged(1:end-1,1);
sum(obs_diffs~=0)
% The values in here are the indices where an observation segment ends.
obs_seg_idx = [find(obs_diffs~=0); num_obs];

obs_seg_idx = [];
for ii = 2:num_obs
    if (ObsMassaged(ii,1) - ObsMassaged(ii-1,1) ~= 0) || ...
            (ObsMassaged(ii,2) - ObsMassaged(ii-1,2) > 60)
        obs_seg_idx = [obs_seg_idx; ii-1];
%         (ObsMassaged(ii,2) - ObsMassaged(ii-1,2))/3600/24
    end
end
obs_seg_idx = [find(ObsMassaged(:,2)/3600/24 <= 50, 1,'last');...
    find(ObsMassaged(:,2)/3600/24 <= 100, 1,'last');...
    find(ObsMassaged(:,2)/3600/24 <= 150, 1,'last');...
    find(ObsMassaged(:,2)/3600/24 <= 200, 1,'last');...
    length(ObsMassaged)];
obs_seg_idx = [12000; length(ObsMassaged)];
% obs_seg_idx = [12000];

% obs_seg_idx = [...
%     find(ObsMassaged(:,2) <= 18674940, 1,'last');...
% %     find(ObsMassaged(:,2) <= 18345000, 1,'last');...
%     find(ObsMassaged(:,2) <= 18863968, 1,'last');...
%     length(ObsMassaged)];

state_len = 6;
IB = [6 6];
state_len = 7;
IB = [7 7];
filter_opts.propagator_opts.OD.A_params.important_block = IB;
filter_opts.important_block = IB;
filter_opts.propagator_opts.OD.state_len = state_len;
filter_opts.SNC_Q = eye(3)*1e-17;
filter_opts.use_SNC = 1;
filter_opts.SNC_Gamma = @(dt) [dt*dt/2 0 0;...
            0 dt*dt/2 0;...
            0 0 dt*dt/2;...
            dt 0 0;...
            0 dt 0;...
            0 0 dt;...
            0 0 0 ];
filter_opts.SRIF_input_P0 = true;
filter_opts.propagator_opts.A_m_ratio =  A_m_ratio; % m2/kg
filter_opts.propagator_opts.A_m_amplitude =  0.02; % m2/kg
filter_opts.propagator_opts.A_m_w =  2*pi/(5*24*3600); % m2/kg
filter_opts.propagator_opts.A_m_amplitude =  0; % m2/kg
filter_opts.propagator_opts.A_m_w =  0; % m2/kg
filter_opts.use_smoother = 0;
P_in = P(1:state_len,1:state_len);


num_seg_iter = 2; % Number of times to iterate the filter on a segment
seg_begin = 1;
iter_state_ap = state_ap(1:state_len);
pfr_total = [];
state_total = [];
cnt = 0;
ell_plot = figure;
% For all segments
for ii = obs_seg_idx'
    cnt = cnt+1;
    if ii ~= seg_begin
        obs_to_process = seg_begin:ii;
    else
        seg_begin = ii + 1;
        continue;
    end
    
    % iterate a few times on the segment
    for jj = 1:num_seg_iter
        [~,X] = ode45(@flyby_two_body_state_dot, ...
            ObsMassaged(obs_to_process,2), ...
            [iter_state_ap; reshape(eye(state_len),IB(1)*IB(2),1)], ...
            filter_opts.ode_opts, filter_opts.propagator_opts);

        filter_opts.ref_state = X;
        myP_in = P_in;

%         if cnt == 2
%             filter_opts.propagator_opts.A_m_ratio = 0.02;
%         else
%             filter_opts.propagator_opts.A_m_ratio =  A_m_ratio;
%         end
        SRIFoutput = SRIF(iter_state_ap, myP_in, ObsMassaged(obs_to_process,:), ...
            filter_opts);
        
%         figure; plot(output.pfr_store(1,:))
        STM_accum = reshape(filter_opts.ref_state(end,state_len+1:end),...
            filter_opts.important_block(1), filter_opts.important_block(2));
        x0_est = STM_accum\SRIFoutput.x_est_store(:,end);
        iter_state_ap = X(1,1:state_len)'+x0_est;
        if ii ~= obs_seg_idx(1)
%             break
        end
    end
    
    if ii == obs_seg_idx(1)
        best_initial_SRIF_state = iter_state_ap;
    end
    
    pfr_total = [pfr_total SRIFoutput.pfr_store];
    state_total = [state_total SRIFoutput.state_store];
    ObsMassaged(seg_begin,1)
    SRIFoutput.state_store(4:6,end)
    SRIFoutput.state_store(7,end)
    find(obs_seg_idx == ii)
    
    output_state = SRIFoutput.state_store(:,end);
    final_P = SRIFoutput.final_P;
    if ii == num_obs
        ii
        color_num = cnt;
        BPlaneTarget_HW7;
        break;
    else
%         if cnt == 1
%             filter_opts.propagator_opts.A_m_ratio = 0.02;
%         else
%             filter_opts.propagator_opts.A_m_ratio =  A_m_ratio;
%         end
        % propagate to next state
        [~,X] = ode45(@flyby_two_body_state_dot, ...
            [ObsMassaged(ii,2) ObsMassaged(ii+1,2)], ...
            [SRIFoutput.state_store(:,end); reshape(eye(state_len),IB(1)*IB(2),1)], ...
            filter_opts.ode_opts, filter_opts.propagator_opts);
        % Propagate the last state
        STM = reshape(X(end,state_len+1:end),...
            filter_opts.important_block(1), ...
            filter_opts.important_block(2));
        iter_state_ap = X(end,1:state_len)';
%         P_in = STM*output.final_P*STM';
        % Keep the covariance positive definite.
        R_bar = SRIFoutput.final_Rj/STM;
        P_in = R_bar;
        filter_opts.SRIF_input_P0 = false;
%         Rinv = inv(R_bar);
%         P_in = Rinv*Rinv';
    end
    seg_begin = ii+1; %Next idx is the next segment
    
    color_num = cnt;
    BPlaneTarget_HW7;
end

% figure(ell_plot);
% legend(plot_handles1, '200-day obs', 'Full obs set');
% xlabel('T (km)')
% ylabel('R (km)')
% title('3\sigma covariance ellipse on B-plane')
% set(gca,'YDir','reverse');

if filter_opts.use_SNC
    plot_title = sprintf('SRIF w/ SNC, Q=%.1e', filter_opts.SNC_Q(1,1));
else
    plot_title = 'SRIF w/o SNC';
end
residual_plot(pfr_total, [0.005, 0.5*1e-6], plot_title)
residual_plot(pfr_total, [0.005, 0.5*1e-6], plot_title, ObsMassaged(1:obs_seg_idx(end),2)/3600/24)
% residual_plot(SRIFoutput.smoothed_pfr, [0.005, 0.5*1e-6], plot_title)
