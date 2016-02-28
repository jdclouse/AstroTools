%% John Clouse IMD HW 4 Problem 2
% 
%% Initialize
clearvars -except hw_pub function_list

V_inf_in=[-5.19425  5.19424 -5.19425];%(km/s)
V_inf_out=[-8.58481  1.17067 -2.42304];%(km/s)
mu_earth=3.986004415e5; %km3/s2

%% Find the B-plane
S_hat = V_inf_in/norm(V_inf_in);
k_hat = [0 0 1]';

T_hat = cross(S_hat,k_hat)/norm(cross(S_hat,k_hat));
R_hat = cross(S_hat, T_hat);
h_hat = cross(V_inf_in,V_inf_out)/norm(cross(V_inf_in,V_inf_out));
B_hat = cross(S_hat, h_hat);

%% The flyby parameters
psi = acos(dot(V_inf_in,V_inf_out)/norm(V_inf_in)/norm(V_inf_out));
rp = mu_earth/(norm(V_inf_in)^2)*(1/cos((pi-psi)/2)-1);

b = mu_earth/(norm(V_inf_in)^2)*sqrt((1+(norm(V_inf_in)^2)*rp/mu_earth)^2 ...
    - 1);
theta = atan2(dot(b*B_hat, R_hat),dot(b*B_hat, T_hat));

%% Results
fprintf('rp = %.3f km\n',rp);
fprintf('psi = %.3f deg\n', psi*180/pi);
fprintf('b = %.3f km\n', b);
fprintf('theta = %.3f deg\n',theta*180/pi);