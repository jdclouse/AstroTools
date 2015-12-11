function [a,e,i,RAAN,w,f] = cart2OE( r, v ,mu)
%cart2OE return classical orbital elements from cartesian coords
% Only valid for e < 1
% units in radians
fcnPrintQueue(mfilename('fullpath')) % Add this code to code app 

h = cross(r,v);
n = cross([0;0;1],h);
ecc_vec = ((norm(v)*norm(v)-mu/norm(r))*r - dot(r,v)*v)/mu;

e = norm(ecc_vec);
a = 0;
if e < 1.0
    specific_energy = norm(v)*norm(v)/2-mu/norm(r);
    a = -mu/2/specific_energy;
end

i = acos(h(3)/norm(h));
RAAN = acos(n(1)/norm(n));
if n(2) < 0
    RAAN = 2*pi-RAAN;
end
w = acos(dot(n,ecc_vec)/(norm(n)*norm(ecc_vec)));
if ecc_vec(3) < 0
    w = 2*pi-w;
end
f = acos(dot(ecc_vec,r)/(norm(ecc_vec)*norm(r)));
if dot(r,v) < 0
    f = 2*pi-f;
end