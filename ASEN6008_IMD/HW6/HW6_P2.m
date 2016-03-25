clear

x0 = 1.2;
x_dot0 = 0;
y0 = 0;
y_dot0 = -1.049657509830343;
X = [x0; y0; 0; x_dot0; y_dot0; 0];

mu = 0.012150585609624;
dunit = 384747.962856037;

T = 6.192169331319632;

[~,X_out] = ode45(@CRTBP, [0,T], X, odeset(),mu);

figure
plot3(X_out(:,1), X_out(:,2), X_out(:,3), 'r')
hold on
axis equal
rad_vec = [0:0.1:2*pi, 2*pi];
my_circ = [cos(rad_vec); zeros(1, length(rad_vec)); sin(rad_vec)]';
for ang = rad_vec
    for blah = 1:length(my_circ)
        new_circ(blah,:) = (Euler2DCM('3', ang)*my_circ(blah,:)')';
    end
    earth =  new_circ * 6378.1/dunit;
    moon =  (new_circ * 1737/dunit);
    plot3(earth(:,1), earth(:,2), earth(:,3))
    plot3(moon(:,1) + 1, moon(:,2), moon(:,3), 'k')
end