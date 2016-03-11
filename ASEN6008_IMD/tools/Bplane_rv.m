function [ b, B_hat, B_plane ] = Bplane_rv( r, v, mu )
%Bplane_rv Compute B Plane from r, v, mu wrt target planet
%   Detailed explanation goes here
fcnPrintQueue(mfilename('fullpath')) % Add this code to code app 

r_mag = norm(r);
v_mag = norm(v);

energy = v_mag*v_mag/2 - mu/r_mag;

v_inf_2 = energy*2;

h_hat = cross(r,v)/norm(cross(r,v));

e_vec = ((v_mag*v_mag-mu/r_mag)*r-dot(r,v)*v)/mu;
a = -mu/v_inf_2;
e = norm(e_vec);

b = -a*sqrt(e*e-1);

phi_s = acos(1/e);
k = [0;0;1];

S_hat = e_vec/e*cos(phi_s) ...
    + cross(h_hat,e_vec)/norm(cross(h_hat,e_vec))*sin(phi_s);
T_hat = cross(S_hat,k)/norm(cross(S_hat,k));
R_hat = cross(S_hat, T_hat);

B_plane = [S_hat T_hat R_hat];

B_hat = cross(S_hat,h_hat);
end

