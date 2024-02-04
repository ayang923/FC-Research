function [barycentric_func] = barycentric_compute_func(w, D_interpolant, phi_interpolant, D)
% Returns polynomial function given barycentric weights and function values
% on interpolation mesh
d = length(w);
denom = 1./(repmat(D, 1, d) - D_interpolant');
barycentric_func = (denom * (phi_interpolant.*w)) ./ (denom*w) ;
end