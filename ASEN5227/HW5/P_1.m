y0 = [1.0; 1.2];

A = [-11 9;9 -11];

h_vec = [0.11;0.10;0.09];

t_end = 2;
figure;
for ii = 1:3
    t = 0;
    h = h_vec(ii);

    Y = y0';
    y = y0;
    T = t;
    while t+h <= t_end
        F = A*y;
        y = y + h*F;
        t = t+h;
        Y = [Y;y'];
        T = [T;t];
    end
    subplot(3,1,ii)
    plot(T,Y(:,1))
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
    subplot(3,1,ii)
    hold on
    plot(t_an,X(:,1),'r')
    ylabel('y_1')
    legend(['Euler, h=' num2str(h_vec(ii))],'Analytical')
end
xlabel('t');
subplot(3,1,1);
    