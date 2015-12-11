function hh = plot_v_time_domain(time, signal, plot_title,xlims, ylims,...
    plot_opts_in)

fcnPrintQueue(mfilename('fullpath')) % Add this code to code app 

plot_opts = '';
if nargin == 6
    plot_opts = plot_opts_in;
end

hh = figure;
plot(time,signal,plot_opts)
title(plot_title)
xlabel('Time (s)')
ylabel('Amplitude (v)') 
grid on

if nargin > 3
    xlim(xlims),ylim(ylims)
end