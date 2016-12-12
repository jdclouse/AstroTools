y0 = [1.0; 1.2];

A = [-11 9;9 -11];

h_vec = [0.11;0.10;0.09];

t_end = 2;
forward_fig = figure;
backwrd_fig = figure;
for ii = 1:3
    t = 0;
    h = h_vec(ii);

    Y_f = y0';
    y = y0;
    T_f = t;
    while t+h <= t_end
        F = A*y;
        y = y + h*F;
        t = t+h;
        Y_f = [Y_f;y'];
        T_f = [T_f;t];
    end
    figure(forward_fig);
    subplot(3,1,ii)
    plot(T_f,Y_f(:,1))
    
    % Backward
    
    t = 0;

    Y_b = y0';
    y = y0;
    T_b = t;
    while t+h <= t_end
        y = (eye(2)-h*A)\y;
        t = t+h;
        Y_b = [Y_b;y'];
        T_b = [T_b;t];
    end
    figure(backwrd_fig);
    subplot(3,1,ii)
    plot(T_b,Y_b(:,1))
end
% Analytical:
v1 = [1;1];
v2 = [1;-1];
c1 = 1.1;
c2 = -0.1;
lam1 = -2;
lam2 = -20;
t_an = 0:0.1:2;

X = zeros(length(t_an),2);
for ii =1:length(t_an)
    X(ii,:) = (c1*exp(lam1*t_an(ii))*v1 + c2*exp(lam2*t_an(ii))*v2)';
end
for ii = 1:3
    for jj = forward_fig:backwrd_fig
        figure(jj);
        subplot(3,1,ii)
        hold on
        plot(t_an,X(:,1),'r')
        ylabel('y_1')
        legend(['Euler, h=' num2str(h_vec(ii))],'Analytical')
    end
end
xlabel('t');
subplot(3,1,1);


figure(forward_fig);
subplot(3,1,1); title('Forward Euler Integration')
figure(backwrd_fig);
subplot(3,1,1); title('Backward Euler Integration')
    