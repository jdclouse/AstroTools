%% John Clouse IMD HW 4 Problem 1
% 
%% Initialize
clearvars -except hw_pub function_list

CelestialConstants

V_spacecraft_wrt_sun=[-10.8559 -35.9372]'; %km/s
v_venus = [15.1945 -31.7927]';%km/s
r_venus = [96948447.3751  46106976.1901]'; %km
mu_sun=1.32712440018e11; %km3/s2
mu_Venus=3.257e5; %km3/s2
R_Venus=    6052; %km

specific_energy_pre = norm(V_spacecraft_wrt_sun)^2/2-mu_sun/norm(r_venus);
fprintf('Heliocentric Energy: %3f km^2/s^2\n',specific_energy_pre);

%% b) turn angle
v_inf = V_spacecraft_wrt_sun - v_venus;
v_inf_2 = norm(v_inf)^2;
num_pts = 200;
rp = linspace(00,200000,num_pts);
turn_angle_store = zeros(1,num_pts);
for ii = 1:num_pts
    turn_angle_store(ii) = pi - 2*acos( 1/(1+v_inf_2*rp(ii)/mu_Venus));
end

figure('Position', hw_pub.figPosn);
plot(rp,turn_angle_store*180/pi,'LineWidth',hw_pub.lineWidth);
title('Turning Angle vs. Planetary Distance')
xlabel('Distance From Venus Center, km')
ylabel('Turning Angle, deg')

%% c) heliocentric energy
energy_leading_pass = zeros(1,num_pts);
energy_trailing_pass = zeros(1,num_pts);
for ii = 1:num_pts
    turn_DCM = Euler2DCM('3',turn_angle_store(ii));
    energy_leading_pass(ii) = norm(turn_DCM(1:2,1:2)*v_inf+v_venus)^2/2 ...
        -mu_sun/norm(r_venus);
    turn_DCM = Euler2DCM('3',-turn_angle_store(ii));
    energy_trailing_pass(ii) = norm(turn_DCM(1:2,1:2)*v_inf+v_venus)^2/2 ...
        -mu_sun/norm(r_venus);
end

figure('Position', hw_pub.figPosn);
hold on
plot(rp,energy_leading_pass,'b','LineWidth',hw_pub.lineWidth);
plot(rp,energy_trailing_pass,'r','LineWidth',hw_pub.lineWidth);
plot([rp(1) rp(end)], ...
    [specific_energy_pre specific_energy_pre],'g-.',...
    'LineWidth',hw_pub.lineWidth);
plot([R_Venus R_Venus], ...
    [max(energy_trailing_pass) min(energy_leading_pass)],'k--',...
    'LineWidth',hw_pub.lineWidth);
legend('Trailing Venus', 'Leading Venus', 'Initial Energy', 'R_{Venus}')
title('Heliocentric Energy After Venus Flyby')
xlabel('Distance From Venus Center, km')
ylabel('Specific Energy, km^2/s^2')

%% d) What does the plot of energy vs. flyby closest approach tell you?
% The spacecraft cannot exit the solar system from just this flyby, as it
% would have to go below the radius of Venus to do so. In addition, the
% energy for both leading and trailing the planet asymptotically approach
% the original energy. This means that the furthest approaches don't add
% much to the heliocentric energy, so trajectories that need to raise their
% aphelian should fly closer.