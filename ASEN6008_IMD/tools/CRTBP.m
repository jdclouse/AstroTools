function state_dot = CRTBP(t, state, mu)

x = state(1);
y = state(2);
z = state(3);
xd = state(4);
yd = state(5);
zd = state(6);
r1 = sqrt((x+mu)^2 + y*y + z*z);
r2 = sqrt((x+mu-1)^2 + y*y + z*z);
ax = -(1-mu)*(x+mu)/r1/r1/r1 - mu*(x-1+mu)/r2/r2/r2 +x +2*yd;
ay = -(1-mu)*y/r1/r1/r1 - mu*y/r2/r2/r2 + y - 2*xd;
az = -(1-mu)*z/r1/r1/r1 - mu*z/r2/r2/r2;

state_dot = [xd;yd;zd;ax;ay;az];