%% HW3 Problem 14
% John Clouse
clear
close all
addpath('C:\Users\John\Documents\Astro\ASEN5050\tools');
%% Initialize

prop_opts.mu = 1;
a_init = 1;
e_init = 0.01;
i_init = pi/6; %30 degrees
RAAN_init = 0;
w_init = pi/4;
f_init = 0;
n_init = sqrt(prop_opts.mu/a_init^3);
[r0, v0 ] = OE2cart( a_init,e_init,i_init,RAAN_init,w_init,f_init,...
    prop_opts.mu);

as = 0.001;

%% Loop over accelerations
P = 2*pi*sqrt(norm(a_init)^3/prop_opts.mu);
ode_opts = odeset('RelTol', 3e-14, 'AbsTol', 1e-20);
X0 = [r0;v0;0;0;as];
% Numerically integrate
[T,X] = ode45(@two_body,[0 P*10], X0, ode_opts, prop_opts);

% get the OEs
[rows, cols] = size(X);
OE = zeros(rows,6);
for jj = 1:rows
    [a,e,i,W,w,f] = cart2OE(X(jj,1:3), X(jj,4:6), prop_opts.mu);
    % make w, W continuous
    if w > pi
        w = w - 2*pi;
    end
    if W > pi
        W = W - 2*pi;
    end
    OE(jj,:) = [a,e,i,W,w,f];
end

% Store the OEs.
OE_store = OE;
T_store = T;
X_store = X;
for ii = 1:length(OE_store(:,6))
    M_store(ii) = E2M(f2E(OE_store(ii,6),OE_store(ii,2)),OE_store(ii,2));
end

%% Plotting
problem='P15';
create_fig_data;

for jj = 1:length(plots)
    figure(plots{jj}.fig);
    if jj <= 2
        subplot(2,1,1)
    end
    plot(T_store,OE_store(:,jj),'LineWidth',hw_pub.lineWidth);
    title(plots{jj}.title,'fontsize',hw_pub.fontSize)
    if jj > 2
        xlabel(plots{jj}.xlabel,'fontsize',hw_pub.fontSize)
    end
    ylabel(plots{jj}.ylabel,'fontsize',hw_pub.fontSize)
end

% analytical delta a plot
E = linspace(0,2*pi*10,1000*10);
E = zeros(1,length(OE_store(:,6)));
for ii = 1:length(OE_store(:,6))
    E(ii) = f2E(OE_store(ii,6),OE_store(ii,2));
end
delta_a = 2*sin(i_init)/n_init^2*(sqrt(1-e_init^2)*cos(w_init)*sin(E)-sin(w_init)*(1-cos(E)))*as;
figure(sma_plots.fig)
subplot(2,1,1)
hold on
plot(T_store,delta_a+a_init,'r--','LineWidth',hw_pub.lineWidth)
legend('Numerical Solution', 'Analytic Solution')
% error plot
subplot(2,1,2)
plot(T_store, abs(OE_store(:,1)-(delta_a+a_init)')./OE_store(:,1)*100,'LineWidth',hw_pub.lineWidth);
ylabel('Percent Error','fontsize',hw_pub.fontSize)
xlabel(plots{1}.xlabel,'fontsize',hw_pub.fontSize)

% analytical delta e plot
for ii = 2:length(OE_store(:,6))
    if E(ii) < E(ii-1)
        E(ii:end) = E(ii:end) + 2*pi;
    end
end
delta_e = sqrt(1-e_init^2)*sin(i_init)/n_init^2/a_init*(3/2*cos(w_init)*E ...
    - 2*e_init*cos(w_init)*sin(E) ...
    + cos(w_init)/4*sin(2*E) ...
    -sqrt(1-e_init^2)/4*sin(w_init)*(1-cos(2*E)))*as;
figure(ecc_plots.fig)
subplot(2,1,1)
hold on 
plot(T_store,delta_e+e_init,'r--','LineWidth',hw_pub.lineWidth)
legend('Numerical Solution', 'Analytic Solution')
% error plot
subplot(2,1,2)
plot(T_store, abs(OE_store(:,2)-(delta_e+e_init)')./OE_store(:,2)*100,'LineWidth',hw_pub.lineWidth);
ylabel('Percent Error','fontsize',hw_pub.fontSize)
xlabel(plots{2}.xlabel,'fontsize',hw_pub.fontSize)

figure;
plot3(X_store(:,1),X_store(:,2),X_store(:,3));
title(sprintf('a_s = %.3f',as));

%% Analytical vs simulated average SMA, eccentricity.
avg_delta_a = 2*sin(i_init)*as/n_init^2*-sin(w_init);
figure('Position', hw_pub.figPosn);
plot(T_store,OE_store(:,1),'LineWidth',hw_pub.lineWidth);
hold on;
plot([T_store(1) T_store(end)],[1 1]*(a_init+avg_delta_a),'r','LineWidth',hw_pub.lineWidth);
% sim_avg_a = sum(OE_store(:,1) ...
%         .*OE_store(:,6))...
%         /sum(OE_store(:,6));
    sim_avg_a = sum(OE_store(:,1))/length(OE_store(:,1));
% sim_avg_a = sum(OE_store(:,1) ...
%         .*M_store')...
%         /sum(M_store);
plot([T_store(1) T_store(end)],[1 1]*(sim_avg_a),'k--','LineWidth',hw_pub.lineWidth);
legend('a','analytic avg \Deltaa', 'numeric avg \Deltaa')
saveas(gcf, ['Figures\' 'P15_avg_a'],'epsc')

% e has secular change. Need to iterate over each orbit to get its avg.
figure('Position', hw_pub.figPosn);
plot(T_store,OE_store(:,2),'LineWidth',hw_pub.lineWidth);
hold on;

X0 = [r0;v0;0;0;as];
T0 = 0;
for num_orbits = 1:10
[T,X] = ode45(@two_body,[T0 P*num_orbits], X0, ode_opts, prop_opts);
    [a,e,i,W,w,f] = cart2OE(X(1,1:3), X(1,4:6), prop_opts.mu);
    avg_delta_e = (3*pi-pi)*as*sqrt(1-e^2)*sin(i)*cos(w)/n_init^2/a;
    
    avg_delta_e = (3/2*pi*cos(w) -sqrt(1-e^2)/4*sin(w)*pi -1/2*cos(w)-sqrt(1-e^2)/2*sin(w))*as*sqrt(1-e^2)*sin(i)/n_init^2/a;
    avg_delta_e = (0)*as*sqrt(1-e^2)*sin(i)/n_init^2/a;
    
    plot([T0 T(end)],[1 1]*(e+avg_delta_e),'r','LineWidth',hw_pub.lineWidth);
    % The true-anomally-weighted avg of eccentricity
    idx = (T_store >= T0).*(T_store < P*num_orbits) == 1;
%     idx = (T_store < P*num_orbits);
    sim_avg_e = sum(OE_store(idx,2) ...
        .*OE_store(idx,6))...
        /sum(OE_store(idx,6));
    sim_avg_e = sum(OE_store(idx,2))/length(OE_store(idx,2));
    plot([T0 T(end)],[1 1]*(sim_avg_e),'k--','LineWidth',hw_pub.lineWidth);
    X0 = X(end,:)';
    T0 = T(end);
    perc_error = (sim_avg_e-(e+avg_delta_e))/(e+avg_delta_e)*100
end

legend('e','analytic avg \Deltae', 'numeric avg \Deltae')
saveas(gcf, ['Figures\' 'P15_avg_e'],'epsc')

% sim_avg_e = sum(OE_store(T_store < 2*pi,2))/length( OE_store(T_store < 2*pi,2))

%% Save!
for jj = 1:length(plots)
    saveas(plots{jj}.fig, ['Figures\' plots{jj}.file],'epsc')
end
