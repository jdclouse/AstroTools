function Idelay = codecarrier_iono_divergence( pseudorange, phase_prange)
%codecarrier_iono_divergence Divergence of pseudorange and phase for
%Ionospheric delay computation
%   Detailed explanation goes here
f_L1 = 1575.42; % MHz
f_L2 = 1227.60; % MHz
c = 2.99792458e8; %m/s
wl_L1 = c/(f_L1*1e6); %m
wl_L2 = c/(f_L2*1e6); %m

Idelay = (pseudorange-phase_prange*wl_L1)/2;
Idelay = Idelay - Idelay(1);
if Idelay(Idelay < 0)
    Idelay = Idelay - min(Idelay);
end

end

