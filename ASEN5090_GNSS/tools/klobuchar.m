function [ Idelay, Iz, t_local ] = klobuchar( az, el, gdlat, gdlon, t, a, b )
%KLOBUCHAR Klobuchar Iono model
%   From Klobuchar, 1987
c = 2.99792458e8; %m/s
A1 = 5e-9; %nighttime value of zenith delay;
A3 = 50400; % s, peak of cosine function, 14h local time
E = el/pi;
A = az/pi;
P=gdlat/pi;
L=gdlon/pi;
t_gps = t; % t is TOW, need TOD

eca = 0.0137./(E+0.11) - 0.022;
[r,col] = size(eca);
if col > r
    eca = eca';
end
PI = P+eca.*cos(A*pi);
if length(PI) < 1
PI = min(0.416, PI);
PI = max(-0.416, PI);
else
    PI(PI > 0.416) = 0.416;
    PI(PI < -0.416) = -0.416;
end

% subionospheric longitude
LI = L + eca.*sin(A*pi)./cos(PI*pi);

% Geomagnetic  lat
Pm = PI + 0.064*cos((LI - 1.617)*pi);

t_local = 4.32e4*LI + t_gps;

OF = 1 + 16*(0.53-E).*(0.53-E).*(0.53-E);

% Pm = Pm*pi;
A2 = a(1) + a(2)*Pm + a(3)*Pm.*Pm + a(4)*Pm.*Pm.*Pm;
if length(A2) < 1
A2 = max([A2 0]);
else
%     if A2(A2 < 0)
%         PI(A2<0)*180
%         A2(A2<0)
%     end
    A2(A2 < 0) = 0;
end
A4 = b(1) + b(2)*Pm + b(3)*Pm.*Pm + b(4)*Pm.*Pm.*Pm;

if length(t_local) > 1
    ind = abs(t_local-A3) < A4/4;
    Iz(ind) = c*(A1 + A2(ind).*cos(2*pi*(t_local(ind)-A3)./A4(ind)));
    Iz(~ind) = A1*c;
    ind = abs(t_local+3600*24-A3) < A4/4;
    Iz(ind) = c*(A1 + A2(ind).*cos(2*pi*(t_local(ind)+3600*24-A3)./A4(ind)));
else
    if abs(t_local-A3) < A4/4 
        Iz = c*(A1 + A2.*cos(2*pi*(t_local-A3)/A4));
    elseif abs(t_local+3600*24-A3) < A4/4
        Iz = c*(A1 + A2.*cos(2*pi*(t_local+3600*24-A3)/A4)); 
    else
        Iz = A1*c;
    end
end

Iz = Iz';
Idelay = Iz.*OF;

end

