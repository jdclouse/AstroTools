gam1 = sym('gam1');
gam2 = sym('gam2');
gam3 = sym('gam3');
gamma_vec = [gam1;gam2;gam3]
gamma_skew = vecSkew(gamma_vec)
omega_body1 = sym('omega_body1');
omega_body2 = sym('omega_body2');
omega_body3 = sym('omega_body3');
omega_bodyma_vec = [omega_body1;omega_body2;omega_body3]
phi = sym('phi');
I = [1,0,0;0,1,0;0,0,1];

B = I+gamma_skew/2 + 1/phi^2*(1-phi/2*cot(phi/2))*gamma_skew*gamma_skew

inv(B)

B_inv = I - (1-cos(phi))/phi^2*gamma_skew + (phi-sin(phi))/phi^3*gamma_skew*gamma_skew

inv(B)== B_inv