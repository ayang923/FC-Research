function [FS_func] = FFT_compute_func(C_n)
% Returns function that computes fourier series with dependency on x mesh
% column vector
%N = (length(C_n) + 1) / 2;
N = length(C_n);
if mod(N, 2) == 0
    k = [0:floor(N/2), -floor(N/2)+1:-1];
else
    k = [0:floor(N/2), -floor(N/2):-1];
end
FS_func = @(x) real(exp(2i*pi*x*k) * C_n);
end

