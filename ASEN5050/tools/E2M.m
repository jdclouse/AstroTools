function M = E2M( E, e )
%E2M Eccentric anomaly (E) to mean anomaly (M)
% units in radians
fcnPrintQueue(mfilename('fullpath')) % Add this code to code app 

% Vallado eqn 2-4
M = E-e*sin(E);