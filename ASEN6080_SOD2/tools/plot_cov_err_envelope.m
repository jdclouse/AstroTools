function plot_handle = plot_cov_err_envelope(cov_store, err_store)

sigma = 2;
len = max(size(cov_store));
cov = arrayfun(@(idx) norm(sqrt(cov_store(:,idx)))*sigma, 1:len);
[rows, num_obs] = size(err_store);
if rows > 1
diff = arrayfun(@(idx) norm(err_store(1:3,idx)), 1:len);
else
    diff = err_store;
end

plot_handle = figure; 
plot(diff);
hold on
plot(cov,'r')
plot(-cov,'r')
title('RSS State error, with covariance envelope')
xlabel('Observation')
ylabel('(m)')
legend('RSS Position Error', '3\sigma covariance envelope')

if rows > 1
RMS =  sqrt(sum(err_store.*err_store,2)./len);
labels = {'x', 'y', 'z'};
figure
for ii = 1:3
subplot(3,1,ii)
hold on
plot(err_store(ii,:),'.')
ylabel(['\delta' labels{ii} '(m)'])
plot(sqrt(cov_store(ii,:))*sigma,'r')
plot(-sqrt(cov_store(ii,:))*sigma,'r')
legend(sprintf('RMS = %f m',RMS(ii)), [num2str(sigma) '\sigma envelope'])

end
xlabel('Observation')
end

if rows == 6
    RMS =  sqrt(sum(err_store.*err_store,2)./len);
    labels = {'vx', 'vy', 'vz'};
    figure
    for ii = 1:3
        subplot(3,1,ii)
        hold on
        plot(err_store(ii+3,:),'.')
        ylabel(['\delta' labels{ii} '(m/s)'])
        plot(sqrt(cov_store(ii+3,:))*sigma,'r')
        plot(-sqrt(cov_store(ii+3,:))*sigma,'r')
        legend(sprintf('RMS = %f m/s',RMS(ii+3)), [num2str(sigma) '\sigma envelope'])

    end
    xlabel('Observation') 
end

perc_in_env = sum(...
    err_store < sqrt(cov_store)*sigma ...
    & err_store > -sqrt(cov_store)*sigma,2)/num_obs*100

perc_outside_env = 100 - perc_in_env