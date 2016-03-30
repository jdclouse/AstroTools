%% John Clouse IMD HW5 problem 1
%% Initialize
clearvars -except hw_pub function_list

AU = 1.49597870691e11; %km
GMs = 1.32712440018e20;
GMem = 4.035032351966808e14;

mu = GMem/(GMs+GMem);

x=15e7;
y = 6e3;
z = 1450;

vx = 0.00075;
vy = 0.08;
vz = 0.019;

t = 450;

period = 2*pi*sqrt(AU*AU*AU/GMs);

%% Results
fprintf('Position:\n')
fprintf('x = %.5e\n', x/(AU/1e3));
fprintf('y = %.5e\n', y/(AU/1e3));
fprintf('z = %.5e\n', z/(AU/1e3));
fprintf('\nVelocity:\n')
fprintf('Vx = %.5e\n', vx/(AU/1e3)*period/2/pi);
fprintf('Vy = %.5e\n', vy/(AU/1e3)*period/2/pi);
fprintf('Vz = %.5e\n', vz/(AU/1e3)*period/2/pi);
fprintf('\nTime = %.4f\n',t/(period/3600/24)*2*pi);

