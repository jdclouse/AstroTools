function state_dot = CRTBP_Halo_Target(t, state, mu)

state_dot = CRTBP(t, state, mu);

x = state(1);
y = state(2);
z = state(3);
xd = state(4);
yd = state(5);
zd = state(6);
r1 = sqrt((x+mu)^2 + y*y + z*z);
r2 = sqrt((x+mu-1)^2 + y*y + z*z);

r1_3 = r1*r1*r1;
r2_3 = r2*r2*r2;
r1_5 = r1_3*r1*r1;
r2_5 = r2_3*r2*r2;

A41 = -(1-mu)/r1_3 + (1-mu)*(x+mu)*3/r1_5*(x+mu) ...
    -mu/r2_3 + mu*(x-1+mu)*3/r2_5*(x+mu-1) + 1;
A42 = (1-mu)*(x+mu)*3/r1_5*y + mu*(x-1+mu)*3/r2_5*y;
A43 = (1-mu)*(x+mu)*3/r1_5*z + mu*(x-1+mu)*3/r2_5*z;

A51 = (1-mu)*y*3/r1_5*(x+mu) + mu*y*3/r2_5*(x+mu-1);
A52 = -(1-mu)/r1_3 + (1-mu)*y*3/r1_5*y -mu/r2_3 + mu*y*3/r2_5*y + 1;
A53 = (1-mu)*y*3/r1_5*z + mu*y*3/r2_5*z;

A61 = (1-mu)*z*3/r1_5*(x+mu) + mu*z*3/r2_5*(x+mu-1);
A62 = (1-mu)*z*3/r1_5*y + mu*z*3/r2_5*y;
A63 = -(1-mu)/r1_3 - mu/r2_3 + (1-mu)*z*z*3/r1_5 + mu*z*z*3/r2_5;

A = [0 0 0 1 0 0;
    0 0 0 0 1 0;
    0 0 0 0 0 1;
    A41 A42 A43 0 2 0;
    A51 A52 A53 -2 0 0;
    A61 A62 A63 0 0 0];

Phi_dot = A*reshape(state(7:end),6,6);

state_dot = [state_dot; reshape(Phi_dot,36,1)];


end
