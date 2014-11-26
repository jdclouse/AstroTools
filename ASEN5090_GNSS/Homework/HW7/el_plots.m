function el_plots(prn_i_want, elevations, times, tstring, sstring, msize)

rows = msize(1);
cols = msize(2);
masked_el = nan(rows,cols);
for ii = 1:length(prn_i_want)
    masked_el(:,prn_i_want(ii)) = elevations(:,prn_i_want(ii));
end
masked_el_vec     = reshape(masked_el,rows*cols,1);

time_mat = repmat(times,1,cols);
time_vec = reshape(time_mat,rows*cols,1);

fig4 = figure; ax4 = axes;
% plot(ax4,time_vec,el_vec,'ob','markerfacecolor','b','markersize',4);
plot(ax4,time_vec,masked_el_vec,'ob','markerfacecolor','b','markersize',4);

ylabel('Elevation (deg)');
xlabel('Time (hr)');
grid(ax4,'on');
when = sprintf('\\bfDuring %s', tstring);
where = sprintf('\\bfat %s', sstring);
title(ax4,{'\fontsize{11}\bfElevation Angle';'\bfof Satellites Seen by Antenna';when;where});