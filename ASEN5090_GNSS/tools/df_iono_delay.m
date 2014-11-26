function [ Idelay, Iz ] = df_iono_delay( rho_L1, rho_L2, el )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
f_L1 = 1575.42; % MHz
f_L2 = 1227.60; % MHz

Idelay = (f_L2*f_L2/(f_L1*f_L1-f_L2*f_L2))*(rho_L2-rho_L1); % Misra/Enge
%If either pseudorange measurement is nonexistent, delay should be a NaN
Idelay(rho_L2 == 0 | rho_L1 == 0) = NaN;
OF = iono_obliq_factor(el);
Iz = Idelay./OF;

end

