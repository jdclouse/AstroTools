%% Problem 4: Constants of motion from Duffing equation
fprintf('Problem 3: Rigid Body Attitude Propagation\n');
clearvars -except function_list pub_opt I I1 I2 I3 T H_sq w_body_store t_mat
close all

A = zeros(3,6);
B = zeros(3,6);
K = zeros(3,6);
w_dot = zeros(3,length(w_body_store),6);

for simulation = 1:6
    fprintf('Simulation %d constants:\n',simulation);
    A(1,simulation) = ((I1-I2)*(2*I3*T(simulation) - H_sq(simulation)) ...
        + (I1-I3)*(2*I2*T(simulation) - H_sq(simulation)))...
        /(I1*I2*I3);
    A(2,simulation) = ((I2-I3)*(2*I1*T(simulation) - H_sq(simulation)) ...
        + (I2-I1)*(2*I3*T(simulation) - H_sq(simulation)))...
        /(I1*I2*I3);
    A(3,simulation) = ((I3-I1)*(2*I2*T(simulation) - H_sq(simulation)) ...
        + (I3-I2)*(2*I1*T(simulation) - H_sq(simulation)))...
        /(I1*I2*I3);
    B(1,simulation) = 2*(I1-I2)*(I1-I3)/I2/I3;
    B(2,simulation) = 2*(I2-I1)*(I2-I3)/I1/I3;
    B(3,simulation) = 2*(I3-I1)*(I3-I2)/I1/I2;
    K(1,simulation) = (2*I2*T(simulation) - H_sq(simulation))...
        *(H_sq(simulation) - 2*I3*T(simulation))/(I1*I1*I2*I3);
    K(2,simulation) = (2*I3*T(simulation) - H_sq(simulation))...
        *(H_sq(simulation) - 2*I1*T(simulation))/(I1*I2*I2*I3);
    K(3,simulation) = (2*I1*T(simulation) - H_sq(simulation))...
        *(H_sq(simulation) - 2*I2*T(simulation))/(I1*I2*I3*I3);
    fprintf('A\n')
    printVector(A(:,simulation),'')
    fprintf('B\n')
    printVector(B(:,simulation),'')
    fprintf('K\n')
    printVector(K(:,simulation),'')
    if simulation <= 3
        fprintf(['With a pure spin about the principle axis, there\n'...
            ,'is no angular acceleration, and it''s like the \n',...
            'spring system starts from zero displacement\n'])
    elseif simulation == 4
        fprintf(['Will remain stable about the 1 axis\n'])
    elseif simulation == 5
        fprintf(['Begins within the linear restoration regime, so cubic\n', ...
            'destabilizing force does not\n'])
    elseif simulation == 6
        fprintf(['Will remain stable about the 3 axis\n'])
    end
    fprintf('\n')
    for idx = 1:length(w_dot)
        w_dot(:,idx,simulation) = ...
            dBodyRatesRigid(w_body_store(:,idx,simulation), I, [0;0;0]);
    end
    figure
    plot(t_mat, w_dot(:,:,simulation).*w_dot(:,:,simulation) ...
        + bsxfun(@times,A(:,simulation),w_body_store(:,:,simulation))...
                                      .*w_body_store(:,:,simulation)...
        + bsxfun(@minus,(bsxfun(@times,B(:,simulation),...
                                w_body_store(:,:,simulation))...
                                .*w_body_store(:,:,simulation)...
                                .*w_body_store(:,:,simulation)...
                                .*w_body_store(:,:,simulation)/2), ...
                                K(:,simulation)))
    plot_title = strcat(['$$\dot{\omega}$$ + $A_{i}\omega_i^2 + B_i\omega_i^4/2 - K_i$ for Simulation ', num2str(simulation)]);
    title(plot_title,'interpreter','latex')
    xlabel('time(s)')
    ylabel('Const. Diffs')
    legend('i = 1', 'i = 2', 'i = 3')
    grid on
end
% + B_i\omega_i^4/2 - K_i

