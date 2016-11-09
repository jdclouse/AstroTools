%% HW3 Problem 13
% John Clouse
clear
close all
addpath('C:\Users\John\Documents\Astro\ASEN5050\tools');
%% Initialize
prop_opts.mu = 1;
r0 = [1;0;0];
v0 = [0;1;0];

accels = [0.001;0.01;0.1;1];

%% Loop over accelerations
P = 2*pi*sqrt(norm(r0)^3/prop_opts.mu);
ode_opts = odeset('RelTol', 3e-14, 'AbsTol', 1e-20);
for ii = 1:length(accels)
    X0 = [r0;v0;accels(ii);0;0];
    % Numerically integrate
    [T,X] = ode45(@two_body,[0 P*10], X0, ode_opts, prop_opts);
    
    % get the OEs
    [rows, cols] = size(X);
    OE = zeros(rows,6);
    for jj = 1:rows
        [a,e,i,w,W,f] = cart2OE(X(jj,1:3), X(jj,4:6), prop_opts.mu);
        OE(jj,:) = [a,e,i,w,W,f];
    end
    
    % Store the OEs in a cell array.
    OE_store{ii} = OE;
    T_store{ii} = T;
    X_store{ii} = X;
end

%% Plotting
problem='P14';
create_fig_data;
f_plots.title    = 'Longitude of periapsis';
plots = {sma_plots, ecc_plots, inc_plots, RAAN_plots, w_plots, f_plots};

for ii = 1:length(accels)
    for jj = 1:length(plots)
        figure(plots{jj}.fig);
        subplot(4,1,ii);
        plot(T_store{ii},OE_store{ii}(:,jj));
        if ii == 1
            title(plots{jj}.title,'fontsize',hw_pub.fontSize)
        end
        if ii == 4
            xlabel(plots{jj}.xlabel,'fontsize',hw_pub.fontSize)
            % zoom for this plot
            figure(sma_plots.fig);
            subplot(4,1,4);
            ylim([-5 5]);
        end
        ylabel(plots{jj}.ylabel,'fontsize',hw_pub.fontSize)
        
    end
end
for jj = 1:length(plots)
    saveas(plots{jj}.fig, ['Figures\' plots{jj}.file],'epsc')
end

r_plots.fig = figure('Position', hw_pub.figPosn);
v_plots.fig = figure('Position', hw_pub.figPosn);
for ii = 1:length(accels)
    % pos
    figure(r_plots.fig);
    subplot(4,1,ii);
    plot(T_store{ii},X_store{ii}(:,1));
    hold on
    plot(T_store{ii},X_store{ii}(:,2),'r');
    if ii == 1
        legend('r_1', 'r_2')
    end
    if ii == 1
        title('Position, DU','fontsize',hw_pub.fontSize)
    end
    if ii == 4
        xlabel('Time, time units','fontsize',hw_pub.fontSize)
    end
    ylabel(sprintf('a_s = %.3f',accels(ii)),'fontsize',hw_pub.fontSize)
    
    % vel
    figure(v_plots.fig);
    subplot(4,1,ii);
    plot(T_store{ii},X_store{ii}(:,4));
    hold on
    plot(T_store{ii},X_store{ii}(:,5),'r');
    if ii == 1
        legend('v_1', 'v_2')
    end
    if ii == 1
        title('Velocity, DU/TU','fontsize',hw_pub.fontSize)
    end
    if ii == 4
        xlabel('Time, time units','fontsize',hw_pub.fontSize)
    end
    ylabel(sprintf('a_s = %.3f',accels(ii)),'fontsize',hw_pub.fontSize)
end
saveas(r_plots.fig, ['Figures\' 'P14_R'],'epsc')
saveas(v_plots.fig, ['Figures\' 'P14_V'],'epsc')

% figure
% plot(T,OE(:,1));
for ii = 1:length(accels)
    figure;
    plot(X_store{ii}(:,1),X_store{ii}(:,2));
    title(sprintf('a_s = %.3f',accels(ii)));
end