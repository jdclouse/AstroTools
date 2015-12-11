function E = M2E( M, e )
%M2E Mean anom (M) to eccentric anom (E)
fcnPrintQueue(mfilename('fullpath')) % Add this code to code app 

tol = 1e-5;

if (M < 0 && M > -pi) || M > pi
    E_1 = M - e;
else
    E_1 = M + e;
end

E = E_1 + tol + 1;
while abs(E_1-E) > tol
    E = E_1;
    E_1 = E - (E - e*sin(E) - M)/(1 - e*cos(E));
end