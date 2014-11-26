%Symbolically Compute Htilda matrix for project
fcnPrintQueue(mfilename('fullpath'))
syms x y z xdot ydot zdot mu J2 Cd s1x s1y s1z s2x s2y s2z s3x s3y s3z
syms Re area m rho theta_dot H r0 rho0 xs ys zs theta
% syms t x(t) y(t) z(t) theta(t)
% xdot = diff(x(t),t);
% ydot = diff(y(t),t);
% zdot = diff(z(t),t);
r = sqrt(x^2+y^2+z^2);
v = sqrt(xdot^2+ydot^2+zdot^2);

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
range = sqrt((x-(xs*cos(theta)-ys*sin(theta)))^2 + ...
    (y-(xs*sin(theta) + ys*cos(theta)))^2 + (z-zs)^2);
% range = sqrt((x(t)-(xs*cos(theta(t))-ys*sin(theta(t))))^2 + ...
%     (y(t)-(xs*sin(theta(t)) + ys*cos(theta(t))))^2 + (z(t)-zs)^2);
range_rate = (x*xdot+y*ydot+z*zdot - (xdot*xs+ydot*ys)*cos(theta) + theta_dot*(x*xs + y*ys)*sin(theta) + (xdot*ys-ydot*xs)*sin(theta) + theta_dot*(x*ys-y*xs)*cos(theta) - zdot*zs)/range;
% range_rate = diff(range,t)
G = [range
    range_rate];
len = 18;
H_tilda = sym('H_tilda',[2 len]);
is_sitex = @(var) var == s1x || var == s2x || var == s3x;
is_sitey = @(var) var == s1y || var == s2y || var == s3y;
is_sitez = @(var) var == s1z || var == s2z || var == s3z;
for ii = 1:2
    for jj = 1:len
        if is_sitex(state(jj))
            H_tilda(ii,jj) = diff(G(ii),xs);
        elseif is_sitey(state(jj))
            H_tilda(ii,jj) = diff(G(ii),ys);
        elseif is_sitez(state(jj))
            H_tilda(ii,jj) = diff(G(ii),zs);
        else
            H_tilda(ii,jj) = diff(G(ii),state(jj));
        end
    end
end