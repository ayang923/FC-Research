function [barycentric_func] = barycentric_compute_func(w, phi_interpolant, D_interpolant)
% Returns polynomial function given barycentric weights and function values
% on interpolation mesh
d = length(w);
denom = @(x) w'./(repmat(x, 1, d) - D_interpolant');
barycentric_func = @(x) sum(denom(x).*phi_interpolant', 2)./sum(denom(x), 2);
end