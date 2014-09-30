%% HW0 Problem 2: Harmonic Oscillator (Numerical Integration)
fprintf('\n');
clearvars -except function_list pub_opt
close all

ode_opts = odeset('RelTol', 1e-12, 'AbsTol', 1e-20);
A = 1.34;
phase = pi/3;
km_rat = 1;
x = [A*cos(phase); -A*sqrt(km_rat)*sin(phase)];

times = linspace(0, 20, 101);

[T,X] = ode45(@harmoscillator, times, x, ode_opts, km_rat);

error = X(:,1)' - analytic_harmonic_oscilator(A, phase, km_rat, times);
plot(times, error);
title('Error between analytical and \newlinenumerically-integrated solns',...
    'fontsize', 16);
xlabel('t', 'fontsize', 14);
ylabel('\Deltax(t)', 'fontsize', 14);
set(gca(),'fontsize',10);
grid on

fprintf('Why would there be an error?\n')
fprintf(...
    ['Ans: The analytical solution is exact and not dependent on \n',...
    '\t previous calculations, while numerical integration is subject\n',...
    '\t to computational error which can compound on previous errors \n',...
    '\t (depending on integrator, time step(s), and model).\n'])