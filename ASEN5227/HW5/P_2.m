close all
clear all

y0 = [1; 0];

alphas = [2,1];
omegas = [1,2];


h_vec = 2/(2+sqrt(3)) + [0.01, 0, -0.01];
h_vec = [h_vec; 1/2 + [0.01, 0, -0.01]];

t_end = 20;
first_fig = figure;
sec_fig = figure;
for ii = 1:length(h_vec)
    for jj = 1:length(alphas)
%     for jj = 1
    t = 0;
    h = h_vec(jj,ii)

%     A = [0 1;-2*alphas(jj) -omegas(jj)^2];
    A = [0 1;-omegas(jj)^2 -2*alphas(jj)];
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
    figure(jj);
    subplot(length(h_vec),1,ii)
    plot(T_f,Y_f(:,1))
    
    % Analitycal
    [V,D] = eigs(A);
    C = V\y0;
    
    v1 = V(:,1);
    v2 = V(:,2);
    c1 = C(1);
    c2 = C(2);
    lam1 = D(1,1);
    lam2 = D(2,2);
    t_an = 0:0.1:t_end;

    X = zeros(length(t_an),2);
    for kk =1:length(t_an)
        X(kk,:) = (c1*exp(lam1*t_an(kk))*v1 + c2*exp(lam2*t_an(kk))*v2)';
    end
    
%     for nn = first_fig:sec_fig
%         figure(nn);
%         subplot(3,1,ii)
        hold on
        plot(t_an,X(:,1),'r')
        ylabel('y_1')
        legend(['Euler, h=' num2str(h_vec(ii))],'Analytical')
%     end
    end
end

%%
% Analytical:
v1 = [1;1];
v2 = [1;-1];
c1 = 1.1;
c2 = -0.1;
lam1 = -2;
lam2 = -20;
t_an = 0:0.1:t_end;

X = zeros(length(t_an),2);
for ii =1:length(t_an)
    X(ii,:) = (c1*exp(lam1*t_an(ii))*v1 + c2*exp(lam2*t_an(ii))*v2)';
end
for ii = 1:3
    for jj = first_fig:backwrd_fig
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
    