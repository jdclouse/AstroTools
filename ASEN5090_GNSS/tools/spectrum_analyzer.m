function hh = spectrum_analyzer(freqs, sig_f, plot_title, xlims, ylims,...
    plot_opts_in)
fcnPrintQueue(mfilename('fullpath')) % Add this code to code app 

plot_opts = '';
if nargin == 6
    plot_opts = plot_opts_in;
end
hh = figure;
plot(freqs,sig_f,plot_opts)
title(plot_title)
xlabel('Frequency (Hz)')
ylabel('Power (dBW)') 
xlim(xlims),ylim(ylims)
grid on