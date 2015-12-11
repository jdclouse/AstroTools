L_boom = 30;
h = L_boom*sin(pi/4)*2;
m_s = 40;
I = m_s*h*h/12
% I = I + ms*2.5^2

m_p = 116; %kg
m = m_p+m_s;
J_s = I; %kg*m2
J_p = 20; %kg*m2
r = 0.88;
s = 0.94;
Bf = 0.79;
Bb = 0.55;
ef = 0.05;
eb=0.55;
rho_s = r*s;
rho_d = (Bf*r*(1-s)+(ef*Bf-eb*Bb)/(ef+eb))*3/2;
P = 4.563e-6; %N/m2
A_sail = h*h;
l = 2; %m
d = l*m_p/m;
b = 0.5; %m

Ft_max = P*A_sail*(1-rho_s);
Fn_max = P*A_sail*(1+rho_s+2/3*rho_d);

a_dd_LHS_1 = J_s + m_s*m_p/m*b*(b+l);
a_dd_LHS_2 = J_p + m_s*m_p/m*l*(b+l);
BIG = a_dd_LHS_2/a_dd_LHS_1;
d_dd_LHS = J_p + m_s*m_p/m*l*(l-b*BIG);

dd_da = m_p/m*(-l+b*BIG)*Ft_max/d_dd_LHS;
dd_dd = -m_p/m*l*Fn_max/d_dd_LHS;

da_da = (-m_p/m*b*Ft_max - m_s*m_p/m*b*l*dd_da)/a_dd_LHS_1;
da_dd = -m_s*m_p/m*b*l*dd_dd/a_dd_LHS_1;

A = [0 1 0 0;
    da_da 0 da_dd 0;
    0 0 0 1;
    dd_da 0 dd_dd 0];

B = [0; -1/a_dd_LHS_1; 0; (1+BIG)/d_dd_LHS];

C = [1 0 0 0];