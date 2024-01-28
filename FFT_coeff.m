function [C_n] = FFT_coeff(phi_D, N)
% Computes coefficients of Fourier Series from discretized grid using FFT
% Assumes phi_D is value of some continuous function on x = 0 to x = 1
% over discrete x mesh.

% phi_D should be a column vector

% Calculates DFT off phi_D, assuming x mesh from x = 0 to x = 1
FFT_phi = fft(phi_D) ./ length(phi_D);

% Trims coefficients for N terms
FFT_phi = FFT_phi(1:N);

% C_n = conjugate(C_{-n})
C_n = [flip(conj(FFT_phi)); FFT_phi(2:end)];
end

