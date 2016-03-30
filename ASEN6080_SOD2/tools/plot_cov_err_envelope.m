function plot_handle = plot_cov_err_envelope(cov_store, err_store)

len = max(size(cov_store));
cov = arrayfun(@(idx) norm(sqrt(cov_store(:,idx)))*3, 1:len);
[rows, ~] = size(err_store);
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
plot(cov_store(ii,:)*2,'r')
plot(-cov_store(ii,:)*2,'r')
legend(sprintf('RMS = %f m',RMS(ii)), '3\sigma envelope')

end
xlabel('Observation')
end