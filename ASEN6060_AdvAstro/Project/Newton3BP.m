%% Newton 3BP
function state_dot = Newton3BP(~, state, opts)

% Assume that mass 3 is negligible. 
% State is: Primary body state, secondary body state, third body state
G = opts.G;
M1 = opts.M1;
M2 = opts.M2;

r1 = state(1:3);
r_dot1 = state(4:6);
r2 = state(7:9);
r_dot2 = state(10:12);
r3 = state(13:15);
r_dot3 = state(16:18);

r_pri_sec = r2-r1;
r_pri_third = r3-r1;
r_sec_third = r3-r2;

% Primary body:
r_dd1 = G*M2/norm(r_pri_sec)^3*r_pri_sec;
r_dd2 = G*M1/norm(r_pri_sec)^3*-r_pri_sec; %Opposite of r_dd1!
r_dd3 = G*M1/norm(r_pri_third)^3*-r_pri_third ...
    + G*M2/norm(r_sec_third)^3*-r_sec_third;


state_dot = [r_dot1; r_dd1; r_dot2; r_dd2; r_dot3; r_dd3];