function [dft,freq] = get_DFT(y,f,Fs,window,overlap)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
freq = 0:Fs/window:Fs/2;
rng = freq >= f(1) & freq <= f(end);
freq = freq(rng);
k = 0:window-1;     % row
n = 0:window-1;     % column
W = exp(-2.*pi.*1i.*n'*k/window);
Windows = buffer(y, window, overlap, 'nodelay');
DFT_mat = W * Windows;
DFT_mat = abs(DFT_mat.^2);
dft = mean(DFT_mat,2);
dft = dft(rng);
end

