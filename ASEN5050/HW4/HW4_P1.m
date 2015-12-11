%% HW4 Problem 1
fprintf('\n');
clearvars -except function_list hw_pub toolsPath 
close all
CelestialConstants; % import useful constants

hp1 = 250; %km
ha1 = 600; %km
hp2 = 2000; %km
ha2 = 5000; %km

a1 = (2*Earth.R + hp1 + ha1)/2;
a2 = (2*Earth.R + hp2 + ha2)/2;

% quick function to compute velocity on the fly:
visviva = @(h,a) sqrt(2*Earth.mu/(Earth.R + h) - Earth.mu/a);
% quick function for pretty output
printout = @(xfer, dv1, dv2, dvtot) fprintf([xfer ' Transfer: ' ...
    'dV1 = ' num2str(dv1, '%.3f') ' km/s, '...
    'dV2 = ' num2str(dv2, '%.3f') ' km/s, '...
    'Total dV = ' num2str(dvtot, '%.3f') ' km/s\n']);

vp1 = visviva(hp1, a1);
va1 = visviva(ha1, a1);
vp2 = visviva(hp2, a2);
va2 = visviva(ha2, a2);
dV_tot_array = [];

% a) initial peri to target apo
a_xfer = (2*Earth.R + hp1 + ha2)/2;
v_xfer_i = abs(visviva(hp1, a_xfer) - vp1);
v_xfer_f = abs(visviva(ha2, a_xfer) - va2);
dV_total = v_xfer_i + v_xfer_f;
printout('P-A', v_xfer_i, v_xfer_f, dV_total);
dV_tot_array = [dV_tot_array dV_total];

% b) initial peri to target peri
a_xfer = (2*Earth.R + hp1 + hp2)/2;
v_xfer_i = abs(visviva(hp1, a_xfer) - vp1);
v_xfer_f = abs(visviva(hp2, a_xfer) - vp2);
dV_total = v_xfer_i + v_xfer_f;
printout('P-P', v_xfer_i, v_xfer_f, dV_total);
dV_tot_array = [dV_tot_array dV_total];

% c) initial apo to target apo
a_xfer = (2*Earth.R + ha1 + ha2)/2;
v_xfer_i = abs(visviva(ha1, a_xfer) - va1);
v_xfer_f = abs(visviva(ha2, a_xfer) - va2);
dV_total = v_xfer_i + v_xfer_f;
printout('A-A', v_xfer_i, v_xfer_f, dV_total);
dV_tot_array = [dV_tot_array dV_total];

% d) initial apo to target peri
a_xfer = (2*Earth.R + ha1 + hp2)/2;
v_xfer_i = abs(visviva(ha1, a_xfer) - va1);
v_xfer_f = abs(visviva(hp2, a_xfer) - vp2);
dV_total = v_xfer_i + v_xfer_f;
printout('A-P', v_xfer_i, v_xfer_f, dV_total);
dV_tot_array = [dV_tot_array dV_total];

[y,i] = min(dV_tot_array);
fprintf(['The minimum dV is for option ' num2str(i) ', total of '...
    num2str(y,'%.3f') ' km/s\n'])