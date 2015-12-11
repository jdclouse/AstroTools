%% HW2 Problem 2
fprintf('\n');
clearvars -except function_list hw_pub toolsPath 
close all

hp = 321; %km
ha = 551; %km
f = 330*pi/180; %rad
delta_t = 65*60; %min
Re = 6378; %km 
mu = 3.986e5; %km3/s2

a = (hp + Re + Re + ha)/2;
e = (a-hp-Re)/a;
n = sqrt(mu/a/a/a);
M0 = E2M(f2E(f,e),e);
M_deploy = M0 + n*delta_t;
if M_deploy > 2*pi
    M_deploy = M_deploy - 2*pi;
end

f_deploy = E2f(M2E(M_deploy,e), e);
fprintf('Shuttle true anom at satellite deployment: %.2f deg\n',...
    f_deploy*180/pi)