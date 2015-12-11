%% HW2 Problem 1
fprintf('\n');
clearvars -except function_list hw_pub toolsPath 
close all
CelestialConstants; % import useful constants

e = 0.15;
a = 2*6378; % km
mu = 3.986e5; % km3/s2
n = sqrt(mu/a/a/a);
d2r = pi/180;

cases = {'A'; 'B'; 'C'; 'D'; 'E'; 'F'};
f = {'0.00'};
E = {'0.00'};
M = {'0.00'};
T = {'0.00'};

f = [f; '30.00'];
E = [E; num2str(f2E(30*d2r,e)/d2r, '%.2f')];
M = [M; num2str(E2M(f2E(30*d2r,e),e)/d2r, '%.2f')];
T = [T; num2str(E2M(f2E(30*d2r,e),e)/n/60, '%.2f')];

f = [f; num2str(E2f(200*d2r,e)/d2r, '%.2f')];
E = [E; '200.00'];
M = [M; num2str(E2M(200*d2r,e)/d2r, '%.2f')];
T = [T; num2str(E2M(200*d2r,e)/n/60, '%.2f')];

f = [f; '90.00'];
E = [E; num2str(f2E(90*d2r,e)/d2r, '%.2f')];
M = [M; num2str(E2M(f2E(90*d2r,e),e)/d2r, '%.2f')];
T = [T; num2str(E2M(f2E(90*d2r,e),e)/n/60, '%.2f')];

f = [f; num2str(E2f(M2E(270*d2r,e),e)/d2r, '%.2f')];
E = [E; num2str(M2E(270*d2r,e)/d2r, '%.2f')];
M = [M; '270.00'];
T = [T; num2str(270*d2r/n/60, '%.2f')];

f = [f; num2str(E2f(E2M(n*25*60,e),e)/d2r, '%.2f')];
E = [E; num2str(E2M(n*25*60,e)/d2r, '%.2f')];
M = [M; num2str(n*25*60/d2r, '%.2f')];
T = [T; '25.00'];


imshow(imread('HW2_P1_graphic1.jpg'));
% disp(table(cases, f, E, M, T))%, 'VariableNames', ...
%     {'Case', 'True Anom (deg)', 'Eccentric Anom (deg)', ...
%     'Mean Anom (deg)', 'T-Tp (min)'}))
% f
% E
% M
% T
