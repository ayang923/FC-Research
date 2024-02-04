function [phi_FT] = FFT_compute_func(C_n, L, D)
% Returns function that computes fourier series with dependency on x mesh
% column vector
%N = (length(C_n) + 1) / 2;
N = length(C_n);
if mod(N, 2) == 0
     freq = [0:N/2, -N/2+1:-1];
else
   freq = [0:(N-1)/2, -(N-1)/2:-1];
end
phi_FT = real(exp(2i*pi*D*freq/L) * C_n);
end

