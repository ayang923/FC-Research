function [w] = barycentric_weights(x_mesh)
% Computes barycentric weights for interpolation dependent only on x_mesh

d = length(x_mesh);

[W_j, W_k] = meshgrid(x_mesh);
W_j(logical(eye(d))) = 1;
W_k(logical(eye(d))) = 0;

w = 1./prod(W_j - W_k, 2);
end

