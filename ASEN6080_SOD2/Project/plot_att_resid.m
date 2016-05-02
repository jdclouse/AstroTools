function plot_att_resid(pfr_store, R, units)
sigs = sqrt(diag(R));
num_obs = length(pfr_store);
figure
offset = 0;
for ii = 1:length(sigs)
    num_subplots = length(sigs);
    if length(sigs) == 6
        num_subplots = 3
    end
    if length(sigs) == 6 && ii == 4
        xlabel('Observation')
        subplot(num_subplots,1,1)
        title('Post-fit Residuals')
        legend('Residual','3\sigma meas envelope')
        figure;
        offset = -3;
    end
    subplot(num_subplots,1,ii+offset)
    plot(1:num_obs, pfr_store(ii,:),'.','LineWidth',1)
    hold on
    plot(1:num_obs,3*sigs(ii)*ones(1,num_obs),'r--')
    plot(1:num_obs,-3*sigs(ii)*ones(1,num_obs),'r--')
%     title(sprintf('Range RMS = %.4e m',output.range_RMS))
    ylabel(units{ii},'Interpreter','latex')
end
xlabel('Observation')
subplot(length(sigs),1,1)
title('Post-fit Residuals')
legend('Residual','3\sigma meas envelope')