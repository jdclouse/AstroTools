t_test = [1:86400]';
a = [3.82e-8
    1.49e-8
    -1.79e-7
    0];
b = [1.43e5
    0
    -3.28e5
    1.13e5];
a = [1.1694e-8 5.9461e-9 -3.6169e-07 4.1277e-07];
b = [1.0511e05 -4.1263e04 -1.3974e06  5.4206e06]';
% 
a = [2.8197e-08  2.4052e-08 -2.7484e-07 -4.7981e-07 ];
b = [1.6588e05 -1.8886e05 -3.3702e06  1.3382e07]';
% for ii = 1:86400
[ Idelay, Iz ] = klobuchar( 210*pi/180*ones(86400,1), 90*pi/180*ones(86400,1), ntus.latgd*pi/180, ntus.lon*pi/180, t_test, a, b );
% [ Idelay, Iz ] = klobuchar( 210*pi/180, 20*pi/180, 40*pi/180, -100*pi/180, t_test, a, b );
figure
plot(t_test/3600,Iz)

el = linspace(0,90)*pi/180;
OF = iono_obliq_factor(el);
figure
plot(el*180/pi, OF)