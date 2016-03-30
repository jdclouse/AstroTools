%% John Clouse IMD HW5 problem 4
%% Initialize
clearvars -except hw_pub function_list
close all


X_ini = [
1.142198291366583
0
-0.1599
0
-0.223
0];

% Constants
mu = 0.012150585609624;
dunit = 384747.962856037;

%% Target a periodic orbit
d = [1;1];
tol = 1e-13;
figure('Position', hw_pub.figPosn)
hold on
while abs(d(1)) > tol && abs(d(2)) > tol
    X = [X_ini; reshape(eye(6),36,1)];

    [T_out,X_out] = ode45(@CRTBP_Halo_Target, [0,2*pi], X, ...
        odeset('Events', @y_crossing),mu);

    d = -[X_out(end,4); X_out(end,6)];
    % STM
    STM = reshape(X_out(end,7:end),6,6);
    y_dot = X_out(end,5);
    state_dot = CRTBP(0,X_out(end,1:6)',mu);

    % The correction while holding x constant
    correction = ([STM(4,3) STM(4,5); STM(6,3) STM(6,5)] ...
        - 1/y_dot*[state_dot(4);state_dot(6)]*[STM(2,3) STM(2,5)])\d;

    X_ini(3) = X_ini(3) + correction(1);
    X_ini(5) = X_ini(5) + correction(2);
    d;
end
[T_out,X_out] = ode45(@CRTBP, [0,T_out(end)*2], X_ini', odeset(),mu);
plot(X_out(:,2), X_out(:,3))
axis equal; xlabel('Y'); ylabel('Z'); title('Targeted periodic orbit, YZ plane')

figure('Position', hw_pub.figPosn)
plot3(X_out(:,1),X_out(:,2),X_out(:,3))
hold on
rad_vec = [0:0.1:2*pi, 2*pi];
my_circ = [cos(rad_vec); zeros(1, length(rad_vec)); sin(rad_vec)]';
for ang = rad_vec
    for blah = 1:length(my_circ)
        new_circ(blah,:) = (Euler2DCM('3', ang)*my_circ(blah,:)')';
    end
    earth =  new_circ * 6378.1/dunit;
    moon =  (new_circ * 1737/dunit);
    plot3(earth(:,1) - mu, earth(:,2), earth(:,3))
    plot3(moon(:,1) + 1-mu, moon(:,2), moon(:,3), 'k')
end
axis equal; xlabel('X'); ylabel('Y'); zlabel('Z'); 
title('Periodic orbit in rotating frame')
fprintf('Initial Conditions:')
fprintf('\tx0 = %.8f\n',X_ini(1));
fprintf('\ty0 = %.8f\n',X_ini(2));
fprintf('\tz0 = %.8f\n',X_ini(3));
fprintf('\tvx0 = %.8f\n',X_ini(4));
fprintf('\tvy0 = %.8f\n',X_ini(5));
fprintf('\tvz0 = %.8f\n',X_ini(6));
