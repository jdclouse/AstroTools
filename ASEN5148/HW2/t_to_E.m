function E = t_to_E(t, e, n)

% t in seconds
% n in rad/s

tol = 1e-6;
E0 = 0;
E1 = 0;
diff = tol + 1;
while abs(diff) > tol
    E0 = E1; % from the last round
    F = (E0-e*sin(E0))/n - t;
    F_prime = (1-e*cos(E0))/n;
    
    E1 = E0 - F/F_prime;
    diff = E1 - E0;
end

E = E1;
% if t > pi/n
%     E = 2*pi-E;
% end