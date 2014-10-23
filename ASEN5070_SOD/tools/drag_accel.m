function accel = drag_accel( state, drag_data )
%drag_accel calculate drag on spacecraft
fcnPrintQueue(mfilename('fullpath')) % Add this code to code app 

Cd = drag_data.Cd;
A = drag_data.A;
m = drag_data.m;

rho0 = 4e-13; %kg/m3
r0 = 7298.145; %km
H = 200.0; %km
theta_dot = 7.29211585530066e-5; %rad/s

if isfield(drag_data, 'model_params')
    if isfield(drag_data.model_params, 'rho0')
        rho0 = drag_data.model_params.rho0;
    end
    if isfield(drag_data.model_params, 'r0')
        r0 = drag_data.model_params.r0;
    end
    if isfield(drag_data.model_params, 'H')
        H = drag_data.model_params.H;
    end
    if isfield(drag_data.model_params, 'theta_dot')
        theta_dot = drag_data.model_params.theta_dot;
    end
end

r = norm(state(1:3));

rho = rho0*exp(-(r-r0)/H);
rel_wind = [state(4) + theta_dot*state(2);
    state(5) - theta_dot*state(1);
    state(6)];

accel = -0.5*Cd*A/m*rho*rel_wind*norm(rel_wind);

end

