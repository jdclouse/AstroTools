function rho = compute_density(r)

rho0 = 4e-13; %kg/m3
r0 = 7298.145; %km
H = 200.0; %km
rho0 = 3.614e-13; %kg/m3
r0 = 700000+6378136.3;
H = 88667;

r_mag = norm(r);

rho = rho0*exp(-(r_mag-r0)/H);