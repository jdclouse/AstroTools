function state = oe2cart(OE)
%cart2oe   Return classical orbital elements given state vector. Earth 
%   orbits only. State given in km, km/s.
%   OE = [a; e; i; RAAN; w; f] = cart2oe(state)
fcnPrintQueue(mfilename('fullpath')) % Add this code to code app 

mu = 3.986e5; % km3/s2

a = OE(1); 
e = OE(2);
i = OE(3);
RAAN = OE(4);
w = OE(5);
f = OE(6);
  
p = a*(1-e*e);
r = p/(1+e*cos(f));
h = sqrt(mu*a*(1-e*e));
heSinf_rp = h*e/r/p*sin(f);

r_vec = zeros(3,1);
r_vec(1) = r*(cos(RAAN)*cos(w+f) - sin(RAAN)*sin(w+f)*cos(i));
r_vec(2) = r*(sin(RAAN)*cos(w+f) + cos(RAAN)*sin(w+f)*cos(i));
r_vec(3) = r*(sin(i)*sin(w+f));

v_vec = zeros(3,1);
v_vec(1) = r_vec(1)*heSinf_rp ...
    - h/r*(cos(RAAN)*sin(w+f) + sin(RAAN)*cos(w+f)*cos(i));
v_vec(2) = r_vec(2)*heSinf_rp ...
    - h/r*(sin(RAAN)*sin(w+f) - cos(RAAN)*cos(w+f)*cos(i));
v_vec(3) = r_vec(3)*heSinf_rp ...
    - h/r*(sin(i)*cos(w+f));

state = [r_vec; v_vec];