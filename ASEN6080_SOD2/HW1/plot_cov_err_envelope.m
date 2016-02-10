function plot_handle = plot_cov_err_envelope(cov_store, err_store)

len = max(size(cov_store));
cov = arrayfun(@(idx) norm(cov_store(:,idx)), 1:len);
diff = arrayfun(@(idx) norm(err_store(1:3,idx)), 1:len);

plot_handle = figure; 
plot(diff);
hold on
plot(cov,'r')
plot(-cov,'r')
title('EKF State error, with covariance envelope')
xlabel('Observation')
ylabel('(m)')
legend('RSS Position Error', '3\sigma covariance envelope')
