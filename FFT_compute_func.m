function [FS_func] = FFT_compute_func(C_n)
% Returns function that computes fourier series with dependency on x mesh
% column vector
N = (length(C_n) + 1) / 2;

FS_func =  @(x) real(exp(2i.*pi.*x*(-N+1:N-1)) * C_n);
end

