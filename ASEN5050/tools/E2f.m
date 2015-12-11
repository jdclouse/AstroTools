function f = E2f( E, e )
%E2f Eccentric anomaly (E) to true anomaly (f)
% Only valid for e < 1
% units in radians
fcnPrintQueue(mfilename('fullpath')) % Add this code to code app 

% Vallado eqn 2-10
f = acos((cos(E) - e)/(1-e*cos(E)));
if E > pi
    f = 2*pi - f;
end