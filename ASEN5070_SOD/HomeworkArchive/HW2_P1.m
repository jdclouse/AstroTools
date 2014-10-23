%% HW2 Problem 1: Orbit propagation with J2 perturbation
%% 1a) Compute the Cartesian partial derivatives of U
% See appendix for hand-derived confirmation of results
fprintf('\n');
clearvars -except function_list pub_opt
close all

syms mu J2 Re x y z
r = (x^2+y^2+z^2)^(1/2);

U = mu/r*(-J2*Re^2/r^2*(3/2*z^2/r^2-1/2));

dU_dx = simplify(diff(U,x))
dU_dy = simplify(diff(U,y))
dU_dz = simplify(diff(U,z))
clearvars -except function_list pub_opt
%% 1b) Propagate orbit, show orbital elements
% All elements but $\Omega$ vary sinusoidally, but the long-term trend 
% is constant (same value at the same point in the period). 
% $\Omega$ decreases with J2 effects.

ode_opts = odeset('RelTol', 1e-12, 'AbsTol', 1e-20);
propagator_opts.J2.use = 1;
r = [-2436.45; -2436.45; 6891.037]; % km
v = [5.088611; -5.088611; 0.0]; % km/s
state = [r;v];

times = 0:20:3600*24;

[T,X] = ode45(@two_body_state_dot, times, state, ode_opts, propagator_opts);

oe_vec = zeros(6,length(times));
for ii = 1:length(times)
    oe_vec(:,ii) = cart2oe(X(ii,:)');
end

plot_oe(oe_vec, times, '')

%% 1c) Compute change in energy

J2 = 0.00108248;
mu = 398600.4; %km3/s2
Re = 6378.145; %km
r_mag = zeros(1,length(times));
v_mag = zeros(1,length(times));
for i = 1:length(times)
    r_mag(i) = norm(X(i,1:3));
    v_mag(i) = norm(X(i,4:6));
end
KE = v_mag.*v_mag/2;%
PE = -3.986e5./r_mag + ...
    mu./r_mag.*(J2*Re*Re./(r_mag.*r_mag).*(3/2*(X(:,3).*X(:,3))'./...
    (r_mag.*r_mag)-1/2));
deltaE = KE + PE - (KE(1)+PE(1));
figure
plot(times/3600, deltaE)
ylabel('\DeltaU (km^2/s^2)')
xlabel('time (hours)')
%% 1d) Compute h_k
h = zeros(length(times),3);
for i = 1:length(times)
   h(i,:) = cross(X(i,1:3), X(i,4:6));
end
figure
plot(times/3600, h(:,3)-h(1,3))
ylabel('\Deltah_z (km^2/s)')
xlabel('time (hours)')