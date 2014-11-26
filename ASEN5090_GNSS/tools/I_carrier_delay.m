function delta_Idelay = I_carrier_delay( phase_rho_L1, phase_rho_L2)
%I_carrier_delay Iono Delay from Carrier differences
%   Detailed explanation goes here
f_L1 = 1575.42; % MHz
f_L2 = 1227.60; % MHz
c = 2.99792458e8; %m/s
wl_L1 = c/(f_L1*1e6); %m
wl_L2 = c/(f_L2*1e6); %m

delta_Idelay = (f_L2*f_L2/(f_L1*f_L1-f_L2*f_L2))*...
    (wl_L1*(phase_rho_L1-phase_rho_L1(1))...
    -wl_L2*(phase_rho_L2-phase_rho_L2(1)));

% Find and correct cycle slips.
% Reset the divergence back to zero when the arbitrary tol is exceeded.
tol = 1; % m
for ii = 2:length(delta_Idelay)
    if abs(delta_Idelay(ii) - delta_Idelay(ii-1)) > tol
        delta_Idelay(ii:end) = delta_Idelay(ii:end) - delta_Idelay(ii);
    end
end
%If either pseudorange measurement is nonexistent, delay should be a NaN
% delta_Idelay(rho_L2 == 0 | rho_L1 == 0) = NaN;
% OF = iono_obliq_factor(el);
% Iz = delta_Idelay./OF;

end

