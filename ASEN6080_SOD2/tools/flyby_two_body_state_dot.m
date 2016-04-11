function state_dot = flyby_two_body_state_dot(t, state, opts)
%two_body_state_dot   Return state_dot given state. Used for numerical
%integration
% The first 6 elements are assumed to be r and v. If a state transition
% matrix (STM) is to be calculated as well, opts.OD needs to be set up.
% Currently assumes r/v are the only state elements that have non-zero
% derivatives.
% fcnPrintQueue(mfilename('fullpath')) % Add this code to code app 

state_dot = zeros(length(state),1);
state_dot(1:3) = state(4:6);
% mu = 3.986e5; % km3/s2
% if isfield(opts, 'mu')
    mu = opts.mu;
% end

r_vec = (state(1:3));
r = norm(r_vec);
state_dot(4:6) = -mu * r_vec/(r*r*r); % Simple 2-body

%% Solar accel
[ r_earth, ~ ] = MeeusEphemeris( opts.Earth, opts.epoch + t/86400 , opts.Sun); %km
r_s_wrt_earth = - opts.EMO2EME*r_earth; %km
delta_s_sc = r_s_wrt_earth - state(1:3); %km
norm_r_s_wrt_earth = norm(r_s_wrt_earth); %km
norm_delta_s_sc = norm(delta_s_sc); %km

state_dot(4:6) = state_dot(4:6) + opts.Sun.mu*...
    (delta_s_sc/(norm_delta_s_sc*norm_delta_s_sc*norm_delta_s_sc) ...
    -r_s_wrt_earth/(norm_r_s_wrt_earth*norm_r_s_wrt_earth*norm_r_s_wrt_earth));

%% SRP
if opts.OD.state_len > 6
srp_param = - opts.solar_flux/opts.c ...
    * opts.au2km*opts.au2km...
    *(opts.A_m_ratio + opts.A_m_amplitude*sin(opts.A_m_w*t))*1e-6;
state_dot(4:6) = state_dot(4:6) ...
    + srp_param...
    *state(7)*delta_s_sc/norm_delta_s_sc/norm_delta_s_sc/norm_delta_s_sc;
else
    srp_param = 0;
end

%% STM
% For orbit determination application
% if isfield(opts, 'OD')
    if opts.OD.use == 1
        
        aug_len = 0;
        
        opts.OD.state_len;
        % The OD.state_len is the length of the estimation state. The rest
        % is the STM, numerically propagated with the A-Matrix
        A = A_state_rvCr(state(1:opts.OD.state_len),opts.OD.A_params, norm_delta_s_sc, r_s_wrt_earth, opts.Sun.mu, srp_param);
        if opts.OD.state_len == 6
            A = A(1:6, 1:6);
        end
        % Block matrix multiplication
%         long_dim = 9;
%         STM = zeros(long_dim);
        STM = reshape(state(opts.OD.state_len+aug_len+1:end),...
            opts.OD.A_params.important_block(1)+0,...
            opts.OD.A_params.important_block(2)+0);
        
        STM_dot = ...
            A(1:opts.OD.A_params.important_block(1),1:opts.OD.A_params.important_block(2))...
            *STM;
        % Pack up the important stuff
        state_dot(opts.OD.state_len+aug_len+1:end) = reshape(STM_dot,...
            (opts.OD.A_params.important_block(1)+0)...
            *(opts.OD.A_params.important_block(2)+0),1);
    end
% end