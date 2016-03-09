function plot_handle = plot_cov_err_envelope2(cov_store, err_store)

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
title('RSS State Deviation, with covariance envelope')
xlabel('Observation')
ylabel('(m)')
legend('RSS Position Deviation', '3\sigma covariance envelope')

if rows > 1
RMS =  sqrt(sum(err_store.*err_store,2)./len);
labels = {'x', 'y', 'z'};
figure
for ii = 1:3
subplot(3,1,ii)
hold on
% plot(err_store(ii,:),'.')
ylabel(['3\sigma ' labels{ii} '(km)'])
plot(sqrt(cov_store(ii,:))*3,'r')
plot(-sqrt(cov_store(ii,:))*3,'r')
% legend(sprintf('RMS = %f m',RMS(ii)), '3\sigma envelope')
end
xlabel('Observation')
subplot(3,1,1)
title('Position Covariance Envelopes')

figure;
for ii = 1:3
subplot(3,1,ii)
hold on
% plot(err_store(ii,:),'.')
ylabel(['3\sigma v' labels{ii} '(km)'])
plot(sqrt(cov_store(ii+3,:))*3,'r')
plot(-sqrt(cov_store(ii+3,:))*3,'r')
% legend(sprintf('RMS = %f m',RMS(ii)), '3\sigma envelope')
end
xlabel('Observation')
subplot(3,1,1)
title('Velocity Covariance Envelopes')


figure;
% plot(err_store(ii,:),'.')
plot(sqrt(cov_store(7,:))*3,'r')
hold on
plot(-sqrt(cov_store(7,:))*3,'r')
% legend(sprintf('RMS = %f m',RMS(ii)), '3\sigma envelope')

ylabel(['3\sigma C_r'])
title('C_r Covariance Envelope')
xlabel('Observation')

end