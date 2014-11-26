%Symbolically Compute A matrix for project
fcnPrintQueue(mfilename('fullpath'))
syms x y z xdot ydot zdot mu J2 Cd s1x s1y s1z s2x s2y s2z s3x s3y s3z
syms Re area m rho theta_dot H r0 rho0
r = sqrt(x^2+y^2+z^2);
v = sqrt(xdot^2+ydot^2+zdot^2);
rel_wind = [xdot + theta_dot*y;
ydot - theta_dot*x;
zdot];
rel_wind_mag = sqrt(rel_wind(1)^2 + rel_wind(2)^2 + rel_wind(3)^2);
state = [x
    y
    z
    xdot 
    ydot 
    zdot 
    mu 
    J2 
    Cd 
    s1x 
    s1y 
    s1z 
    s2x 
    s2y 
    s2z 
    s3x 
    s3y 
    s3z];
F = [xdot
    ydot
    zdot
    (-mu/r^3*x + 1.5*mu*J2*Re^2/(r^5)*(5*z^2/r^2 - 1)*x -0.5*Cd*area/m*rho0*exp(-(r-r0)/H)*rel_wind_mag*rel_wind(1))
    (-mu/r^3*y + 1.5*mu*J2*Re^2/(r^5)*(5*z^2/r^2 - 1)*y -0.5*Cd*area/m*rho0*exp(-(r-r0)/H)*rel_wind_mag*rel_wind(2))
    (-mu/r^3*z + 1.5*mu*J2*Re^2/(r^5)*(5*z^2/r^2 - 3)*z -0.5*Cd*area/m*rho0*exp(-(r-r0)/H)*rel_wind_mag*rel_wind(3))
    0
    0
    0
    0
    0
    0
    0
    0
    0
    0
    0
    0];
len = 18;
A = sym('A',[len len]);
for ii = 1:len
    for jj = 1:len
        A(ii,jj) = diff(F(ii),state(jj));
    end
end

% diff(state, F)
