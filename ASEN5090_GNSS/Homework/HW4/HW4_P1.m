%% HW4: GPS Positioning Accuracy
%% Initialize
fprintf('\n');
clearvars -except function_list pub_opt
close all

%% 1) Ellipsoidal height vs time
% The variation of ellipsoidal height seems to vary more when fewer
% satellites are used in the position solution.
load('GPS2009L1L2_data.mat');
time = data(:,1);
lat = data(:,2);
long = data(:,3);
ell_h = data(:,4);
num_sat = data(:,5);
[haxes, hline1, hline2] = plotyy(time, ell_h, time, num_sat);
ylim(haxes(2), [min(num_sat)-1, max(num_sat)+1]);
xlabel('GPS Time of Week (s)');
ylabel(haxes(1), 'Ellipsoidal height (m)');
ylabel(haxes(2), 'Number of Satellites');

%% 2) Latitude/Longitude errors
% The scatter looks more varied in the y direction (North/South) than the x
% direction.
mean_lat = mean(lat);
lat_err = lat - mean_lat;
mean_long = mean(long);
long_err = long - mean_long;
nn = length(lat_err);
sig_lat = std(lat);
sig_mean_lat = sig_lat/sqrt(nn);
sig_long = std(long);
sig_mean_long = sig_long/sqrt(nn);
fprintf(sprintf('Mean latitude: %%.%df degrees\n', ...
    dec_places(sig_mean_lat)),mean_lat);
fprintf(sprintf('Mean longitude: %%.%df degrees\n', ...
    dec_places(sig_mean_long)),mean_long);

% Get the error arc-lengths
re = 6378e3; % m, assuming spherical
NS_err = 2*pi*re*lat_err/360; 
EW_err = cos(lat*pi/180)*2*pi*re.*long_err/360; 

figure
plot(EW_err, NS_err, '.')
title('Position Errors')
ylabel('North (m)')
xlabel('East (m)')
grid on
axis equal
hold on

%% 3) Compute standard deviations, P, error ellipse, 50% CEP

% Standard deviations:
sig_ell_h = std(ell_h);
fprintf('Latitude standard deviation = %.0e m\n', std(NS_err));
fprintf('Longitude standard deviation = %.0e m\n', std(EW_err));
fprintf('Ellipsoidal height standard deviation = %.0e m\n', sig_ell_h);

P = cov([EW_err NS_err]);
fprintf('Covariance matrix:\n');
for i = 1:2
    fprintf('\t')
    for j = 1:2
        if P(i,j) > 1
            fprintf('%4.3f\t', P(i,j));
        else
            fprintf(sprintf('%%.%df\t', dec_places(P(i,j))+3),P(i,j));
        end
    end
    fprintf('\n')
end

% Eigenvalues of the covariance matrix are the principle components
[evec,ev]=eig(P);
ell_a=sqrt(ev(1,1));%, evec(:,1)'
ell_b=sqrt(ev(2,2));%, evec(:,2)'
angle=atan2(evec(2,1),evec(1,1));
drawellipse(ell_a,ell_b,angle,0,0,'g-');
fprintf('Error Ellipse semimajor axis = %.3e m\n', ell_a);
fprintf('Error Ellipse semiminor axis = %.3e m\n', ell_b);

err_radii = sqrt((NS_err.*NS_err)+(EW_err.*EW_err));
CEP_50_radius = err_radii(round(nn/2));
fprintf('CEP 50%% radius = %.3e m\n', CEP_50_radius);
drawellipse(CEP_50_radius,CEP_50_radius,0,0,0,'r-.');