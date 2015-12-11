function E = newtonRaphson( func, deriv, vars, guess, tol )
%newtonRaphson Newton-Raphson method, taking function and derivative
% derivative function handles
fcnPrintQueue(mfilename('fullpath')) % Add this code to code app 

x = guess;
x_1 = x + tol + 1; % ensures at least one iteration
while abs(x_1-x) > tol
    x_1 = x - func(vars)/deriv(vars);
end