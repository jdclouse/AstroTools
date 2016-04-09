
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

state_len = 6;
IB = [6 6];
state_len = 7;
IB = [7 7];
filter_opts.propagator_opts.OD.A_params.important_block = IB;
filter_opts.important_block = IB;
filter_opts.propagator_opts.OD.state_len = state_len;
P_in = P(1:state_len,1:state_len);

num_seg_iter = 4; % Number of times to iterate the filter on a segment
seg_begin = 1;
iter_state_ap = state_ap(1:state_len);
pfr_total = [];
state_total = [];
% For all segments
for ii = obs_seg_idx'
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
        if seg_begin == 1 || (ObsMassaged(seg_begin,2) - ObsMassaged(seg_begin-1,2))/3600/24 < 5
            myP_in = P_in;
        else
            myP_in = P_in*10;
        end
        output = SRIF(iter_state_ap, myP_in, ObsMassaged(obs_to_process,:), ...
            filter_opts);
        
%         figure; plot(output.pfr_store(1,:))
        STM_accum = reshape(filter_opts.ref_state(end,state_len+1:end),...
            filter_opts.important_block(1), filter_opts.important_block(2));
        x0_est = STM_accum\output.x_est_store(:,end);
        iter_state_ap = X(1,1:state_len)'+x0_est;
    end
    
    if ii == obs_seg_idx(1)
        best_initial_SRIF_state = iter_state_ap;
    end
    
    pfr_total = [pfr_total output.pfr_store];
    state_total = [state_total output.state_store];
    ObsMassaged(seg_begin,1)
    output.state_store(4:6,end)
    output.state_store(7,end)
    find(obs_seg_idx == ii)
    if ii == num_obs
        ii
        break;
    else
        % propagate to next state
        [~,X] = ode45(@flyby_two_body_state_dot, ...
            [ObsMassaged(ii,2) ObsMassaged(ii+1,2)], ...
            [output.state_store(:,end); reshape(eye(state_len),IB(1)*IB(2),1)], ...
            filter_opts.ode_opts, filter_opts.propagator_opts);
        % Propagate the last state
        STM = reshape(X(end,state_len+1:end),...
            filter_opts.important_block(1), ...
            filter_opts.important_block(2));
        iter_state_ap = X(end,1:state_len)';
%         P_in = STM*output.final_P*STM';
    end
    seg_begin = ii+1; %Next idx is the next segment
end

residual_plot(pfr_total, [0.005, 0.5*1e-6], '1 iter')
residual_plot(pfr_total, [0.005, 0.5*1e-6], 'SRIF', ObsMassaged(:,2))

% [~,X] = ode45(@flyby_two_body_state_dot, ObsMassaged(obs_to_process,2), ...
%     [state_ap; reshape(eye(7),49,1)], ...
%         filter_opts.ode_opts, filter_opts.propagator_opts);
%     
% filter_opts.ref_state = X;
% output_1 = SRIF(state_ap, P, ObsMassaged(obs_to_process,:), filter_opts);
% % figure
% % plot3(X(:,1),X(:,2),X(:,3))
% % xlabel('x (km)')
% % ylabel('y (km)')
% % zlabel('z (km)')
% % axis equal
% 
% % The last point at which the residuals looked good.
% iter_point = 231; %5300
% if length(ObsMassaged(obs_to_process,2)) < iter_point
%     iter_point = length(ObsMassaged(obs_to_process,2));
% end