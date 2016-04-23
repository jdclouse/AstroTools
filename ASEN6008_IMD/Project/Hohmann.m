%%Hohmann
v_earth_circ = sqrt(Sun.mu/Earth.a);
v_jup_circ = sqrt(Sun.mu/Jupiter.a);

a_xfer = (Earth.a + Jupiter.a)/2;
vp_xfer = sqrt(2*Sun.mu/Earth.a - Sun.mu/a_xfer);
va_xfer = sqrt(2*Sun.mu/Jupiter.a - Sun.mu/a_xfer);

T = pi*sqrt(a_xfer^3/Sun.mu)/3600/24/365.25;