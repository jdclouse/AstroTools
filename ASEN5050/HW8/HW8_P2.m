%% HW8 Problem 2
fprintf('\n');
clearvars -except function_list hw_pub toolsPath 
close all
CelestialConstants; % import useful constants

planets = {Mercury, Venus, Earth, Moon, Mars, Jupiter, Saturn,...
    Uranus, Neptune};
line_style = {'g', 'm', 'b', 'k--', 'r', 'r-.', 'b--',...
    'k', 'c'};
h = 800;
figure('Position',[0 0 hw_pub.figWidth hw_pub.figHeight])
for idx = 1:length(planets)
    planet = planets{idx};
    ss_rate{idx} = 360/planet.P_days; %deg/day
    fprintf(['Sun-synchronous nodal regression rate for ' ...
        planet.name ...
        ': %.3e deg/Earth day\n'], ss_rate{idx})
    r = h + planet.R;
    p = r;
    n = sqrt(planet.mu/r^3)*180/pi*day2sec;
    i_vec = 0:0.01:180;
    nodal_regression_rate = -3*n*planet.R^2*planet.J2/2/p^2*cosd(i_vec);
    reg_rate_cell_array{idx} = nodal_regression_rate;
    plot(i_vec, nodal_regression_rate, line_style{idx},'LineWidth',2)
    hold on
    target_i{idx} = acosd(ss_rate{idx}*2*p^2/(-3*n*planet.R^2*planet.J2));
    names{idx} = planet.name;
end
legend(names, 'Location', 'northwest')
xlabel('Inclination (deg)')
ylabel('$\dot{\Omega}$ (deg/Earth day)','interpreter', 'latex')

for idx = 1:length(planets)
    planet = planets{idx};
    fprintf(['\n' planet.name ': Target Inclination: '])
    if ss_rate{idx} > max(reg_rate_cell_array{idx})
        fprintf('Not achievable')
    else
        fprintf('%.2f degrees', target_i{idx})
    end
end
fprintf('\n')