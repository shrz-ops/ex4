function [PowerSpectrum,freq] = get_FFT(y,f,Fs)
% Summary of this function goes here
%   Detailed explanation goes here
L = length(y);      % Length of signal
freq = 0:Fs/L:Fs/2;
rng = freq >= f(1) & freq <= f(end);
freq = freq(rng);
yft = fft(y);
yft_norm = yft/length(yft);
PowerSpectrum = abs(yft_norm) .^ 2;
PowerSpectrum = PowerSpectrum(rng);
end

