function plot_att_resid(pfr_store, R, units)
sigs = sqrt(diag(R));
num_obs = length(pfr_store);
figure
for ii = 1:length(sigs)
    subplot(length(sigs),1,ii)
    plot(1:num_obs, pfr_store(ii,:),'.','LineWidth',1)
    hold on
    plot(1:num_obs,3*sigs(ii)*ones(1,num_obs),'r--')
    plot(1:num_obs,-3*sigs(ii)*ones(1,num_obs),'r--')
%     title(sprintf('Range RMS = %.4e m',output.range_RMS))
    ylabel(units{ii})
end
xlabel('Observation')