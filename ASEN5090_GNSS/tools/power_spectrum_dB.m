function P = power_spectrum_dB(signal)
P = 20*log10(abs(fft(signal)));