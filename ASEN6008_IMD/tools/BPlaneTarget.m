function [b, B_hat, B_plane, psi, rp] = BPlaneTarget(V_inf_in, V_inf_out, mu)
fcnPrintQueue(mfilename('fullpath')) % Add this code to code app 

V_inf_in_2 = norm(V_inf_in)^2;

S_hat = V_inf_in/norm(V_inf_in);
k_hat = [0 0 1]';

T_hat = cross(S_hat,k_hat)/norm(cross(S_hat,k_hat));
R_hat = cross(S_hat, T_hat);
h_hat = cross(V_inf_in,V_inf_out)/norm(cross(V_inf_in,V_inf_out));
B_hat = cross(S_hat, h_hat);

B_plane = [S_hat T_hat R_hat];

% Turning angle
psi = acos(dot(V_inf_in,V_inf_out)/norm(V_inf_in)/norm(V_inf_out));

% Closes approach to planet
rp = mu/V_inf_in_2*(1/cos((pi-psi)/2)-1);

% magnitude of B vector
b = mu/V_inf_in_2*sqrt((1+V_inf_in_2*rp/mu)^2 - 1);