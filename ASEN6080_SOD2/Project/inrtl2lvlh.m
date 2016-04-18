function DCM = inrtl2lvlh(r,v)
% r = state(1:3);
% v = state(4:6);
z = -r/norm(r);
h = cross(r,v);
y = -h/norm(h);
x = cross(y,z)/norm(cross(y,z));

DCM = [x y z]';