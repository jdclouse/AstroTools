clear

IC_set(:,1) = [
    -0.08
-0.03
0.01
3.5
-3.1
-0.1
26];

IC_set(:,2) = [0.05
-0.05
0
4.0
2.6
0
25];

IC_set(:,3) = [0.8300
0
0.114062816271683
0
0.229389507175582
0
15];

IC_set(:,4) = [-0.05
-0.02
0
4.09
-5.27
0
15];

for ii = 1:4
    X = IC_set(1:end-1,ii);

    mu = 0.012150585609624;
    dunit = 384747.962856037;

    T = 6.192169331319632;

    [T_out,X_out] = ode45(@CRTBP, [0,IC_set(end,ii)], X, odeset(),mu);

    figure
    subplot(2,1,1);
    plot3(X_out(:,1), X_out(:,2), X_out(:,3), 'r')
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
    axis equal
    
    X_inrt = X_out;
    X_inrt(:,1) = X_inrt(:,1) + mu;
    for jj = 1:length(T_out)
        t = T_out(jj);
        ang = t;
        X_inrt(jj,1:3) = (Euler2DCM('3', -ang)*X_inrt(jj,1:3)')';
    end
    subplot(2,1,2);
    plot3(X_inrt(:,1), X_inrt(:,2), X_inrt(:,3), 'r')
    hold on
    for ang = rad_vec
    for blah = 1:length(my_circ)
        new_circ(blah,:) = (Euler2DCM('3', ang)*my_circ(blah,:)')';
    end
    earth =  new_circ * 6378.1/dunit;
    plot3(earth(:,1) - mu, earth(:,2), earth(:,3))
    end
    plot3(my_circ(:,1), my_circ(:,3), my_circ(:,2), 'k')
    axis equal
    
end