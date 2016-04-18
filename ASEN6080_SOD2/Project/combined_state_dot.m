function state_dot = combined_state_dot(t, state, opts)

state_dot = zeros(13,1);
state_dot(1:6) = two_body_state_dot(t, state(1:6), opts.PV_prop_opts);

% att_prop_opts.DCM_inrtl2lvlh = inrtl2lvlh(state(1:3), state(4:6))...
%     *state(1:3);
DCM_inrtl2body = EulerParam2DCM(state(7:10));

opts.att_prop_opts.R_body = DCM_inrtl2body*state(1:3);
state_dot(7:end) = rigid_body_quat(t, state(7:end), opts.att_prop_opts);