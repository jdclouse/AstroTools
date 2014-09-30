function E = get_eccentric_anomaly(f, e)
%get_eccentric_anomaly Calculate eccentric anomaly
fcnPrintQueue(mfilename('fullpath')) % Add this code to code app 

E = acos((e + cos(f))./(1+e.*cos(f)));

E(f > pi) = 2*pi - E(f > pi);
