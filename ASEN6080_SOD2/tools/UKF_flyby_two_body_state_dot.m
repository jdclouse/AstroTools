function state_dot = UKF_flyby_two_body_state_dot(t, state, opts)
%UKF_two_body_state_dot Calculate state derivative for UKF sigma points
% fcnPrintQueue(mfilename('fullpath')) % Add this code to code app 

opts.OD.use = 0; % Don't use this

% Propagate each sigma point. 
% The incoming state is 2L+1 separate states of length L
state_dot = zeros(length(state),1);
L = opts.UKF.L;
for ii = 1:(2*L+1)
    idx_1 = L*(ii-1)+1;
    idx_L = idx_1 + L - 1;
    state_dot(idx_1:idx_L) = ...
        flyby_two_body_state_dot(t, state(idx_1:idx_L), opts);
end

end

