function fig = residual_plot(pfr, sigs, the_title)

num_obs = length(pfr);

range_RMS = sqrt(sum(pfr(1,:).*pfr(1,:))/num_obs);
rangerate_RMS = sqrt(sum(pfr(2,:).*pfr(2,:))/num_obs);

fig = figure;
subplot(2,1,1)
plot(pfr(1,:),'.','LineWidth',1)
hold on
plot(3*sigs(1)*ones(1,num_obs),'r--')
plot(-3*sigs(1)*ones(1,num_obs),'r--')
title([the_title ' Range Residuals'])
ylabel('m')
legend(sprintf('RMS = %.4e m',range_RMS),'3\sigma measurement noise')
subplot(2,1,2)
plot(pfr(2,:),'.','LineWidth',1)
hold on
plot(3*sigs(2)*ones(1,num_obs),'r--')
plot(-3*sigs(2)*ones(1,num_obs),'r--')
title([the_title ' Range-Rate, '])
ylabel('m/s'),xlabel('Observation')
legend(sprintf('RMS = %.4e m/s',rangerate_RMS),'3\sigma measurement noise')


subplot(2,1,1)