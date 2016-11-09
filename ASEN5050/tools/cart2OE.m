function [a,e,i,RAAN,w,f] = cart2OE( r, v ,mu)
%cart2OE return classical orbital elements from cartesian coords
% Only valid for e < 1
% units in radians
% fcnPrintQueue(mfilename('fullpath')) % Add this code to code app 

e_zero_tol = 1e-10;
i_zero_tol = 1e-12;
n_zero_tol = 1e-12;

h = cross(r,v);
n = cross([0;0;1],h);
ecc_vec = ((norm(v)*norm(v)-mu/norm(r))*r - dot(r,v)*v)/mu;

e = norm(ecc_vec);
a = 0;
if e < 1.0
    specific_energy = norm(v)*norm(v)/2-mu/norm(r);
    a = -mu/2/specific_energy;
elseif e > 1.0    
    specific_energy = norm(v)*norm(v)/2-mu/norm(r);
    a = -mu/2/specific_energy;
end

i = acos(h(3)/norm(h));

% Check if we are near-equatorial.
if abs(i) < i_zero_tol || norm(n) <  n_zero_tol
    RAAN = 0;
else
    RAAN = acos(n(1)/norm(n));
    if n(2) < 0
        RAAN = 2*pi-RAAN;
    end
end

% Check if we are near-circular
if e < e_zero_tol || abs(i) < i_zero_tol || norm(n) <  n_zero_tol
    w = 0;
else
    w = acos(dot(n,ecc_vec)/(norm(n)*norm(ecc_vec)));
    if ecc_vec(3) < 0
        w = 2*pi-w;
    end
end

% f is always assigned here, even if not the true anomally
% the special cases are always detectable through i, w, and RAAN
if e > e_zero_tol && i > i_zero_tol % inclined, elliptical
    f = acos(dot(ecc_vec,r)/(norm(ecc_vec)*norm(r)));
    if dot(r,v) < 0
        f = 2*pi-f;
    end

    % Special case if f is imaginary
    if imag(f)
        f = r(1)/norm(r); %use arg of latitude
        if r(2)<0
            f = 2*pi-f;
        end
        w = 0;
    end
elseif e < e_zero_tol && i > i_zero_tol % inclined, circular
    % argument of latitude u = w+f
    u = acos(dot(n,r)/norm(n)/norm(r));
    if r(3) < 0
        u = 2*pi-u;
    end
    f = u; % Assign to the true anomaly slot.
elseif e > e_zero_tol && i < i_zero_tol % Equatorial, elliptical
    % longitude of periapsis = w+RAAN
    lop = acos(ecc_vec(1)/e);
    if ecc_vec(2) < 0
        lop = 2*pi-lop;
    end
    f = lop; % Assign to the true anomaly slot.
elseif e < e_zero_tol && i < i_zero_tol % Equatorial, circular
    % true longitude lambda = lop + M
    lambda = acos(r(1)/norm(r));
    if r(2) < 0
        lambda = 2*pi-lambda;
    end
    f = lambda; % Assign to the true anomaly slot.
end
    