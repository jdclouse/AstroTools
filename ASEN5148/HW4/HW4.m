%% HW 4: Thermal
%% Problem 1
% 1. Calculate the worst case hot and cold temperatures for a 1 meter 
% diameter spherical aluminum spacecraft in 600km altitude Earth equatorial
% orbit.
% ?s = ?IR = 0.25 specific heat = 900 J/Kg-°K

% Assume s/c is directly between earth and sun for worst-case hot.
% Assume earth is between sun and s/c for worst-case cold.
% Assume no internal heat generated.

stef_boltz = 5.67e-8; % W/(m2K4)
c = 900;
absorptivity = 0.25;
emissivity = 0.25;

G_sun = 1371; % W/m2 (Brown)
q_earthIR = 237; % W/m2 (Brown)
albedo = 0.3; % ratio of reflected sunlight (Brown)

d = 1; % m, diameter
A_projected = pi*d*d/4;
A_surf = pi*d*d;

Re = 6378; % km
H = 600; % km
F = 0.5*(1-sqrt(H*H + 2*H*Re)/(H+Re));
Ka = 0.657 + 0.54*Re/(Re+H) - 0.196*(Re/(Re+H))^2;

Q_sun = G_sun*A_projected*absorptivity;
Q_albedo = G_sun*albedo*A_surf*absorptivity*F*Ka ;
Q_IR = q_earthIR*emissivity*A_surf*F; % abs. = emiss. at a given wvlngth
Q_internal = 0;

T_hot = ((Q_sun + Q_albedo + Q_IR+Q_internal)...
    /(stef_boltz*emissivity*A_surf))^(1/4)
T_cold = ((Q_IR+Q_internal)/(stef_boltz*emissivity*A_surf))^(1/4)

% Calculate the same for a 1 m x 1 m x 25mm solar array honeycomb panel 
% with solar cells on one side, and white paint on the other.

% Assume solar cells face sun for worst-case hot
% Assume panel full-on faces the sun/earth
% Assume paint faces earth for worst-case cold

clearvars -except stef_boltz c G_sun q_earthIR albedo Re H Ka

solar_cell_abs = 0.805;
solar_cell_emi = 0.825;
white_paint_abs = 0.252;
white_paint_emi = 0.853;

A_projected = 1; % m2
A_surf = 1+1+4*0.025; %m2

F = Re*Re/(H+Re)^2;

Q_sun = G_sun*A_projected*solar_cell_abs;
Q_albedo = G_sun*albedo*A_projected*white_paint_abs*F*Ka; %paint facing earth
Q_IR = q_earthIR*white_paint_emi*A_projected*F; %paint facing earth

T_hot = ((Q_sun + Q_albedo + Q_IR)/...
    (stef_boltz*(solar_cell_emi+white_paint_emi)*A_projected))^(1/4)
T_cold = ((Q_IR)/(stef_boltz*(solar_cell_emi+white_paint_emi)*A_projected))^(1/4)

%% Problem 2
% a. Calculate the average temperature of the spacecraft at 5 minute 
% orbital intervals for the first 5 orbits. Assume a 100kg spherical 
% spacecraft made of aluminum launched into a sun-synchronous Earth orbit 
% with a starting temperature of 12°C, with 300 watts internal dissipation.
% b. Plot the family of curves for 2.a. for 0.05 increments of ?IR from 
% 0.25 to 0.9
clearvars -except stef_boltz c G_sun q_earthIR albedo Re H Ka
close all

absorptivity = 0.25;
emissivity = 0.25;

d = 1; % m, diameter
A_projected = pi*d*d/4;
A_surf = pi*d*d;

C2K = 273; % K, Celsius to Kelvin
Q_internal = 300; % W
m = 100; % kg

mu = 3.986e5; % km3/s2
n = sqrt(mu/(Re+H)^3);
P = 2*pi*sqrt((Re+H)^3/mu);
F = 0.5*(1-sqrt(H*H + 2*H*Re)/(H+Re));

alb_frac = 0.5; % Assuming dawn/dusk sun-sync orbit

Q_sun = G_sun*A_projected*absorptivity;
Q_albedo = G_sun*albedo*alb_frac*A_surf*absorptivity*F*Ka;
Q_IR = q_earthIR*emissivity*A_surf*F; % abs. = emiss. at a given wvlngth

eclipse_entry = pi/2 - asin(Re/(Re+H));
eclipse_exit = pi/2 + asin(Re/(Re+H));

dt = 5*60;
duration = P*5;
for emissivity = 0.25:0.05:0.9
    theta = -pi/2; %starting at the equator. 
    % assuming it's going north on the sun side
    T0 = 12 + C2K; % k
    T_array = zeros(length(0:dt:duration),1);
    ii = 1;
    for t = 0:dt:duration
        if theta >= eclipse_entry && theta <= eclipse_exit
            Q_sun = 0;
            Q_albedo = 0;
        else
            alb_frac = 1;
            Q_sun = G_sun*A_projected*absorptivity;
            Q_albedo = G_sun*albedo*alb_frac*A_surf*absorptivity*F*Ka;
        end
        dT_dt = (Q_sun + Q_albedo + Q_IR + Q_internal ...
            - emissivity*A_surf*stef_boltz*T0^4)/m/c;
        T_array(ii) = T0 ;
        T0 = T0 + dT_dt*dt;
        theta = theta + n*dt;
        if theta > pi % keep it between -pi and pi
            theta = theta - 2*pi;
        end
        ii = ii + 1;
    end
    plot((0:dt:duration)/60,T_array-C2K);
    hold on
end
xlabel('Time (minutes)')
ylabel('Spacecraft Temp (degrees C)')
