%% HW7 Problem 1
fprintf('\n');
clearvars -except function_list hw_pub toolsPath 
close all
CelestialConstants; % import useful constants

X0 = [5492.000;%km
 3984.001 ;%km
 2.955 ;%km
 -3.931 ;%km/sec
 5.498 ;%km/sec
 3.665 ];%km/sec

% Anon fcn to calculate specific energy. It shouldn't change!
spec_energy = @(X) norm(X(4:6))^2/2 - Earth.mu/norm(X(1:3));

% Classical orbit elements
[a,e,i,RAAN,w,f] = cart2OE(X0(1:3),X0(4:6),Earth.mu);

% Get the stuff that's propagated
n = sqrt(Earth.mu/a/a/a);
M0 = E2M(f2E(f,e),e);

for t = [100 1e6]; %s
    % Final mean anom is easy...
    Mf = M0 + n*t;
    % Unwinde the mean anom
    while Mf > 2*pi
        Mf = Mf - 2*pi;
    end
    % Final true anom
    ff = E2f(M2E(Mf,e),e);
    % Back to ECI!
    [r_f, v_f] = OE2cart(a,e,i,RAAN,w,ff,Earth.mu);
    fprintf('r_f(t=%d):\n',t)
    disp(r_f);
    fprintf('delta Energy(t=%d):\n',t)
    disp(spec_energy([r_f;v_f]) - spec_energy(X0));
end

% Anonymous function to calculate 2-body accel
two_body = @(t,X) [X(4);X(5);X(6);...
    -Earth.mu*X(1)/norm(X(1:3))^3;...
    -Earth.mu*X(2)/norm(X(1:3))^3;...
    -Earth.mu*X(3)/norm(X(1:3))^3];

% Anon fcn to calculateposition difference.
calc_dr = @(X_exp, r_f) sqrt((X_exp(1)-r_f(1))^2 ...
        +(X_exp(2)-r_f(2))^2 ...
        +(X_exp(3)-r_f(3))^2);

tol=1e-12;
options=odeset('RelTol',tol,'AbsTol',[tol tol tol tol tol tol]);
for t = [100 1e6]
    [t_array,X_array]=ode45(two_body,[0 t],X0,options);
    fprintf('r_f(t=%d)(integrated):\n',t_array(end))
    disp(X_array(end,1:3)');
    fprintf('delta Energy(t=%d)(integrated):\n',t)
    disp(spec_energy(X_array(end,1:6)) - spec_energy(X0));
end

for tol = [1e-12 1e-10 1e-8 1e-6 1e-4]
    options=odeset('RelTol',tol,'AbsTol',[tol tol tol tol tol tol]);
    [t_array,X_array]=ode45(two_body,[0 1e6],X0,options);
    fprintf('Position Diff @ tol = %e:\n',tol)
    disp(calc_dr(X_array(end,1:3)',r_f));
    fprintf('delta Energy @ tol = %f:\n',tol)
    disp(spec_energy(X_array(end,1:6)) - spec_energy(X0));
end
