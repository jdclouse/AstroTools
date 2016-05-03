function plot_att_resid(pfr_store, R, units)

hw_pub.figWidth = 1120; % pixels
hw_pub.figHeight = 840; % pixels
hw_pub.figPosn = [0, 0, hw_pub.figWidth, hw_pub.figHeight];
% Example: some_fig = figure('Position', hw_pub.figPosn);
hw_pub.lineWidth = 2; % pixels
hw_pub.fontSize = 12;

sigs = sqrt(diag(R));
num_obs = length(pfr_store);
figure('Position', hw_pub.figPosn);
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
        figure('Position', hw_pub.figPosn);
        offset = -3;
    end
    subplot(num_subplots,1,ii+offset)
    plot(1:num_obs, pfr_store(ii,:),'.','LineWidth',1)
    hold on
    plot(1:num_obs,3*sigs(ii)*ones(1,num_obs),'r--')
    plot(1:num_obs,-3*sigs(ii)*ones(1,num_obs),'r--')
%     title(sprintf('Range RMS = %.4e m',output.range_RMS))
    ylabel(units{ii},'Interpreter','latex','fontsize',hw_pub.fontSize)
end
xlabel('Observation','fontsize',hw_pub.fontSize)
subplot(num_subplots,1,1)
title('Post-fit Residuals','fontsize',hw_pub.fontSize)
h = legend('Residual','3\sigma meas envelope');
set(h,'fontsize',hw_pub.fontSize);