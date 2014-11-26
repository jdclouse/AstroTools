function HW8_P2_plot(eps, diffs, row, filter, logx_fig, loglog_fig)
figure(logx_fig)
subplot(4,2,(row-1)*2+1)
semilogx(eps, diffs(1,:))
ylabel(sprintf('%s X1', filter))
if row == 4
   xlabel('\epsilon')
end
subplot(4,2,(row-1)*2+2)
semilogx(eps, diffs(2,:))
ylabel(sprintf('%s X2', filter))
if row == 4
   xlabel('\epsilon')
end

figure(loglog_fig)
subplot(4,2,(row-1)*2+1)
loglog(eps, abs(diffs(1,:)))
ylabel(sprintf('%s X1', filter))
if row == 4
   xlabel('\epsilon')
end
subplot(4,2,(row-1)*2+2)
loglog(eps, abs(diffs(2,:)))
ylabel(sprintf('%s X2', filter))
if row == 4
   xlabel('\epsilon')
end


end