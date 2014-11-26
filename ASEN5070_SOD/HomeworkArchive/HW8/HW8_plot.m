function HW8_plot(eps, diffs, subplot_sector, filter, logx_fig, loglog_fig)
figure(logx_fig)
subplot(2,2,subplot_sector)
semilogx(eps, diffs)
title(filter)
xlabel('\epsilon')
ylabel('Trace Difference')

figure(loglog_fig)
subplot(2,2,subplot_sector)
loglog(eps, abs(diffs)) %ISSUE

title(filter)
xlabel('\epsilon')
ylabel('Trace Difference')
end