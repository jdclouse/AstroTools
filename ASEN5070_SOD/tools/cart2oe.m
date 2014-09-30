function OE = cart2oe(state)
%cart2oe   Return classical orbital elements given state vector. Earth 
%   orbits only. State given in km, km/s.
%   OE = [a; e; i; RAAN; w; f] = cart2oe(state)
fcnPrintQueue(mfilename('fullpath')) % Add this code to code app 

mu = 3.986e5; % km3/s2
r_vec = state(1:3);
r = norm(r_vec);
v_vec = state(4:6);
v = norm(v_vec);
h = cross(r_vec,v_vec);

% Specific Energy:
E = v*v/2 - mu/r;

a = -mu/2/E;
e = sqrt(1-dot(h,h)/a/mu);
i = acos(h(3)/norm(h));
RAAN = atan2(h(1), -h(2));

arg_lat = atan2(r_vec(3)/sin(i), (r_vec(1)*cos(RAAN)+r_vec(2)*sin(RAAN)));
cosf = (a*(1-e*e)-r)/e/r;
f = acos(max([min([1, cosf]), -1]));
if dot(r_vec,v_vec) < 0
    f = 2*pi - f;
end
w = arg_lat - f;
if w < 0
    w = w + 2*pi;
end

OE = [a; e; i; RAAN; w; f];