function output = analytic_harmonic_oscilator(amplitude, phase, ...
    ang_freq_2, time_array)
%analytic_harmonic_oscilator   Return output for harmonic oscilator given
%   amplitude, phase, (angular frequency)^2, and array of time.
fcnPrintQueue(mfilename('fullpath')) % Add this code to code appendix

output = amplitude*cos(sqrt(ang_freq_2)*time_array + phase);