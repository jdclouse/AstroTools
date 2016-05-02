
ang_label = {'$\psi$ (deg)', '$\theta$ (deg)', '$\phi$ (deg)'};
rate_label = {'$\dot{\psi}$ (deg/s)', '$\dot{\theta}$ (deg/s)', '$\dot{\phi}$ (deg/s)'};
ang_err_label = {'$\delta\psi$ (deg)', '$\delta\theta$ (deg)', '$\delta\phi$ (deg)'};
rate_err_label = {'$\delta\dot{\psi}$ (deg/s)', '$\delta\dot{\theta}$ (deg/s)', '$\delta\dot{\phi}$ (deg/s)'};

figure; 
for ii = 1:3
subplot(3,1,ii)
plot(filter_output.state_store(ii,:)*180/pi)
hold on
plot(X_out_angles(:,ii+6)'*180/pi, 'r')
ylabel(ang_label{ii},'Interpreter','latex');
end
xlabel('Observation')
subplot(3,1,1);
title(sprintf([case_title ' Angles, SNC = %.1e'],filter_opts.SNC_Q(1)))
legend('Estimate','Truth')

% rates
figure; 
for ii = 1:3
subplot(3,1,ii)
plot(filter_output.state_store(ii+3,:)*180/pi)
hold on
plot(rate_meas_data(:,1+ii)*180/pi, 'r')
ylabel(rate_label{ii},'Interpreter','latex');
end
xlabel('Observation')
subplot(3,1,1);
title(sprintf([case_title ' Rates, SNC = %.1e'],filter_opts.SNC_Q(1)))
legend('Estimate','Truth')

% Error
figure; 
num_pts = length(filter_output.state_store(ii,:));
for ii = 1:3
subplot(3,1,ii)
plot((filter_output.state_store(ii,:)-euler_angs(1:num_pts,ii)')*180/pi)
% plot((filter_output.state_store(ii,:)-X_out_angles(1:num_pts,ii+6)')*180/pi)
hold on
plot(3*sqrt(filter_output.cov_store(ii,:))*180/pi,'r')
plot(-3*sqrt(filter_output.cov_store(ii,:))*180/pi,'r')
ylabel(ang_err_label{ii},'Interpreter','latex');
end
xlabel('Observation')
subplot(3,1,1);
title(sprintf([case_title ' Angle Error, SNC = %.1e'],filter_opts.SNC_Q(1)))
legend('Error','3\sigma envelope')

figure; 
num_pts = length(filter_output.state_store(ii,:));
for ii = 1:3
subplot(3,1,ii)
plot((filter_output.state_store(ii+3,:)-sim_out(1:num_pts,ii+10)')*180/pi)
% plot((filter_output.state_store(ii+3,:)-X_out_angles(1:num_pts,ii+6+3)')*180/pi)
hold on
plot(3*sqrt(filter_output.cov_store(ii+3,:))*180/pi,'r')
plot(-3*sqrt(filter_output.cov_store(ii+3,:))*180/pi,'r')
ylabel(rate_err_label{ii},'Interpreter','latex');
end
xlabel('Observation')
subplot(3,1,1);
title(sprintf([case_title ' Rate Error, SNC = %.1e'],filter_opts.SNC_Q(1)))
legend('Error','3\sigma envelope')

% Residuals
plot_att_resid(filter_output.pfr_store, filter_opts.R, resid_plot_units)