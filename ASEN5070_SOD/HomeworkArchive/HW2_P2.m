%% HW2 Problem 2: Orbit propagation with J2 perturbation and drag
%% 2a) Propagate orbit, compute change in energy
% Energy is not conserved because drag is not a conservative force.
propagator_opts.drag.use = 1;
propagator_opts.drag.Cd = 2.0;
propagator_opts.drag.A = 3.6;
propagator_opts.drag.m = 1350;

[T2,X2] = ode45(@two_body_state_dot, times, state, ode_opts, ...
    propagator_opts);

oe_vec2 = zeros(6,length(times));
for ii = 1:length(times)
    oe_vec2(:,ii) = cart2oe(X2(ii,:)');
end

%Change in energy
r_mag = zeros(1,length(times));
v_mag = zeros(1,length(times));
for i = 1:length(times)
    r_mag(i) = norm(X2(i,1:3));
    v_mag(i) = norm(X2(i,4:6));
end
KE = v_mag.*v_mag/2;
PE = -3.986e5./r_mag ...
    + mu./r_mag.*(J2*Re*Re./(r_mag.*r_mag).*(3/2*(X(:,3).*X(:,3))'./...
    (r_mag.*r_mag)-1/2));
deltaE = KE + PE - (KE(1)+PE(1));
figure
plot(times, deltaE)
ylabel('\Deltah_z (km^2/s^2)')
xlabel('time (hours)')

%% 2b) Plot the OE diffs from Problem 1
% The semi-major axis reduces as a result of the drag.
plot_oe(oe_vec, times, '\Delta',...
    2*pi/sqrt(mu/(oe_vec2(1)*oe_vec2(1)*oe_vec2(1))), oe_vec2)