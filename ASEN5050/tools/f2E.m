function E = f2E( f, e )
%f2E True anomaly (f) to eccentric anomaly (E)
% Only valid for e < 1
% units in radians
fcnPrintQueue(mfilename('fullpath')) % Add this code to code app 

% Vallado eqn 2-9
E = acos((e + cos(f))/(1+e*cos(f)));
if f > pi
    E = 2*pi - E;
end