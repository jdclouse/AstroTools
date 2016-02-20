function fig = residual_plot(pfr, sigs)

fig = figure;
subplot(2,1,1)
plot(1:num_obs, pfr(1,:),'.','LineWidth',1)
hold on
plot(1:num_obs,3*sigs(1)*ones(1,num_obs),'r--')
plot(1:num_obs,-3*sigs(1)*ones(1,num_obs),'r--')
% title(sprintf('Range RMS = %.4e m',output.range_RMS))
ylabel('m')
subplot(2,1,2)
plot(1:num_obs, pfr(2,:),'.','LineWidth',1)
hold on
plot(1:num_obs,3*sigs(2)*ones(1,num_obs),'r--')
plot(1:num_obs,-3*sigs(2)*ones(1,num_obs),'r--')
% title(sprintf('Range-Rate RMS = %.4e m/s',output.rangerate_RMS))
ylabel('m/s'),xlabel('Observation')


subplot(2,1,1)