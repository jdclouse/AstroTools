function OF = iono_obliq_factor( el )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
Re = 6378.1370e3;
h_I = 350e3;
zenith_ang = pi/2-el;

OF = 1./sqrt(1-(Re*sin(zenith_ang)/(Re+h_I)).*(Re*sin(zenith_ang)/(Re+h_I)));
end

