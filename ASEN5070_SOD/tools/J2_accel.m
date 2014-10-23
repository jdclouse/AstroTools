function accel = J2_accel( pos, params )
%J2_accel Acceleration due to J2
fcnPrintQueue(mfilename('fullpath')) % Add this code to code app 

% Defaults (earth)
J2 = 0.00108248;
mu = 398600.4; %km3/s2
Re = 6378.145; %km

% Use input params if desired
if nargin > 1 % There are params
    if isfield(params, 'J2')
        J2 = params.J2;
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
const = 1.5*mu*J2*Re*Re/(r*r*r*r*r);
sin_sq_phi = z*z/(r*r);

accel = const*[5*sin_sq_phi - 1;
    5*sin_sq_phi - 1;
    5*sin_sq_phi - 3].*pos;
end

