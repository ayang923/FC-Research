function [C_n] = FFT_coeff(phi_D)
% Computes coefficients of Fourier Series from discretized grid using FFT
% Assumes phi_D is value of some continuous function on x = 0 to x = 1
% over discrete x mesh, not including x=1;

% phi_D should be a column vector

% Calculates DFT off phi_D, assuming x mesh from x = 0 to x = 1
C_n = fft(phi_D) ./ length(phi_D);
end

