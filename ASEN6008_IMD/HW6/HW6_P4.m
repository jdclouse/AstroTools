X_ini = [
1.142198291366583
0
-0.1599
0
-0.223
0];

d = [1;1]
tol = 1e-13;
figure
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

correction = ([STM(4,3) STM(4,5); STM(6,3) STM(6,5)] ...
    - 1/y_dot*[state_dot(4);state_dot(6)]*[STM(2,3) STM(2,5)])\d;

X_ini(3) = X_ini(3) + correction(1);
X_ini(5) = X_ini(5) + correction(2);
d
plot(X_out(:,2), X_out(:,3))
end
axis equal