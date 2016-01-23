function accel = J3_accel( pos, params )
%J3_accel Acceleration due to J3
fcnPrintQueue(mfilename('fullpath')) % Add this code to code app 

% Defaults (earth)
J3 = -2.5327e-6;
mu = 398600.4; %km3/s2
Re = 6378.145; %km

% Use input params if desired
if nargin > 1 % There are params
    if isfield(params, 'J3')
        J3 = params.J3;
    end
    if isfield(params, 'mu')
        mu = params.mu;
    end
    if isfield(params, 'Re')
        Re = params.Re;
    end
end

%Calculate accel
r = norm(pos);
z = pos(3);
const = 2.5*mu*J3*Re*Re*Re/(r*r*r*r*r*r*r);
sin_sq_phi_z = z*z*z/(r*r);

accel = const*[7*sin_sq_phi_z - 3*z;
    7*sin_sq_phi_z - 3*z;
    7*sin_sq_phi_z - 6*z + 0.6*r*r/z].*pos;
end

