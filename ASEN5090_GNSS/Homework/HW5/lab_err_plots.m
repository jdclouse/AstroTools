function lab_err_plots(x,y,x_label, y_label, ptitle, wpts)

plot(x, y)
hold on
ylims = [min(y) max(y)];
hold on
for ii = 1:length(wpts)
    plot([wpts(ii),wpts(ii)],ylims,'r--')
    text(wpts(ii),ylims(2),sprintf('WP%d',ii))
end
hold off
datetick('x', 'HH:MM')
title(ptitle, 'FontWeight', 'bold')
ylabel(y_label)
xlabel(x_label)