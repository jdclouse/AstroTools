function f = fpa_to_f(fpa, e)

% t in seconds
% n in rad/s

tol = 1e-6;
f0 = 0;
f1 = 0;
diff = tol + 1;
while abs(diff) > tol
    f0 = f1; % from the last round
    F = tan(fpa) + e*tan(fpa)*cos(f0) - sin(f0);
    F_prime = -e*tan(fpa)*sin(f0) - cos(f0);
    
    f1 = f0 - F/F_prime;
    diff = f1 - f0;
end

f = f1;