%% HW0 Problem 1: Harmonic Oscillator (Analytic)
fprintf('\n');
clearvars -except function_list pub_opt
close all

times = linspace(0, 20, 101);
plot(times, analytic_harmonic_oscilator(1.34, pi/3, 1, times));
xlabel('t', 'fontsize', 16);
ylabel('x(t)', 'fontsize', 16);
set(gca(),'fontsize',12);
grid on