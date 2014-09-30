%% Problem 3: Rigid Body Attitude Propagation
fprintf('Problem 3: Rigid Body Attitude Propagation\n');
clearvars -except function_list pub_opt
close all
I1 = 125;
I2 = 100;
I3 = 75;
I = [I1 0 0; 0 I2 0; 0 0 I3];

% In Euler Params
Attitude = [1;0;0;0];
body_w = zeros(3,6);
body_w(:,1) = [1;0;0]; % rad/s
body_w(:,2) = [0;1;0]; % rad/s
body_w(:,3) = [0;0;1]; % rad/s
body_w(:,4) = [1;0.1;0]; % rad/s
body_w(:,5) = [0;1;0.1]; % rad/s
body_w(:,6) = [0.1;0;1]; % rad/s

discussion = cell(1,6);
discussion{1} = ['Simulation 1:\n'...
                 'There are no disturbance torques. w1 does not change\n',...
                 'because w2 and w3 are zero.\n'...
                 'Will remain in a pure spin about this axis (minimum\n'...
                 'energy case)\n\n'];
discussion{2} = ['Simulation 2:\n'...
                 'There are no disturbance torques. w2 does not change\n',...
                 'because w3 is zero.\n'...
                 'Will remain in a pure spin about this axis \n'...
                 '(intermediate energy case, no disturbance)\n\n'];
discussion{3} = ['Simulation 3:\n'...
                 'There are no disturbance torques. w3 does not change\n',...
                 'because w2 is zero.\n'...
                 'Will remain in a pure spin about this axis (maximum\n'...
                 'energy case)\n\n'];
discussion{4} = ['Simulation 4:\n'...
                 'There are no disturbance torques. w1 changes as w2\n',...
                 'is non-zero, and as w3 becomes non-zero.\n'...
                 'Spin is fairly stable about the I1 axis.\n\n'];
discussion{5} = ['Simulation 5:\n'...
                 'There are no disturbance torques. w2 changes as w3\n',...
                 'is non-zero, and as w1 becomes non-zero.\n'...
                 'Spin is not stable due to initial conditions located\n'...
                 'along the sepratrix.\n\n'];
discussion{6} = ['Simulation 6:\n'...
                 'There are no disturbance torques. w3 changes as w1\n',...
                 'is non-zero, and as w2 becomes non-zero.\n'...
                 'Spin is fairly stable about the I3 axis.\n\n'];

% Propagation time
end_time = 60; % seconds
delta_t = 0.1; % seconds
t_mat = 0:delta_t:end_time;
[rows, cols] = size(t_mat);
% For use in problem 4:
T = zeros(6,1); %Kinetic energy
H_sq = zeros(6,1); %Angular momentum squared
w_body_store = zeros(3,cols,6);

for simulation = 1:6
    % Arrays for recording and plotting
    t_mat = 0:delta_t:end_time;
    [rows, cols] = size(t_mat);
    EP_mat = zeros(4,cols);
    EP_mat(:,1) = Attitude;
    H = I*body_w(:,simulation);
    T(simulation,1) = 0.5*dot(H,body_w(:,simulation));
    H_sq(simulation,1) = dot(H,H);
    w_body_store(:,1,simulation) = body_w(:,simulation);
    idx = 2;

    % Begin propagation of initial conditions
    state = [Attitude; body_w(:,simulation)];
    for t = 0:delta_t:end_time-delta_t
        % Euler integration
%         delta_state = [dEulerParam(state(1:4), state(5:7)); ...
%             dBodyRatesRigid(state(5:7), I, [0;0;0])] * delta_t;
%         state = state + delta_state;
        
        % RK4 integration
        k1 = [dEulerParam(state(1:4), state(5:7)); ...
            dBodyRatesRigid(state(5:7), I, [0;0;0])];
        k2 = [dEulerParam(state(1:4) + delta_t*k1(1:4)/2, state(5:7) + delta_t*k1(5:7)/2); ...
            dBodyRatesRigid(state(5:7) + delta_t*k1(5:7)/2, I, [0;0;0])];
        k3 = [dEulerParam(state(1:4) + delta_t*k2(1:4)/2, state(5:7) + delta_t*k2(5:7)/2); ...
            dBodyRatesRigid(state(5:7) + delta_t*k2(5:7)/2, I, [0;0;0])];
        k4 = [dEulerParam(state(1:4) + delta_t*k3(1:4), state(5:7) + delta_t*k3(5:7)); ...
            dBodyRatesRigid(state(5:7) + delta_t*k3(5:7), I, [0;0;0])];
        
        state = state + delta_t/6*(k1 + 2*k2 + 2*k3 + k4);

        % Normalize EP
        state(1:4) = state(1:4)/norm(state(1:4));
        
        EP_mat(:,idx) = state(1:4);
        
        
        H = I*state(5:7);
%         T(simulation,idx) = 0.5*dot(H,state(5:7));
%         H_sq(simulation,idx) = dot(H,H);
        w_body_store(:,idx,simulation) = state(5:7);
        
        idx = idx + 1;
    end

    figure
    plot(t_mat, EP_mat);
    plot_title = strcat(['EP Elements for Simulation ', num2str(simulation)]);
    title(plot_title)
    xlabel('time(s)')
    ylabel('Element Magnitude')
    legend('\beta_{0}', '\beta_{1}', '\beta_{2}', '\beta_{3}')
    grid on
    fprintf(discussion{simulation})
%     figure
%     plot(t_mat,H_sq(simulation,:)-H_sq(simulation,1))
%     title('delta H_sq')
%     figure
%     plot(t_mat,T(simulation,:)-T(simulation,1))
%     title('delta T')
end