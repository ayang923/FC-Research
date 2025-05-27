clear; close all;

f = @(x, y) x.^2-y.^2;

%smooth case
l_1 = @(theta) cos(2*pi*theta)+0.65*cos(4*pi*theta)-0.65;
l_2 = @(theta) sin(2*pi*theta);
l_1_prime = @(theta) -2*pi*sin(2*pi*theta) - 0.65*4*pi*sin(4*pi*theta);
l_2_prime = @(theta) 2*pi*cos(2*pi*theta);
l_1_dprime = @(theta) -4*pi^2*cos(2*pi*theta) - 0.65*(4*pi)^2*cos(4*pi*theta);
l_2_dprime = @(theta) -4*pi^2*sin(2*pi*theta);

K_boundary = @(theta_1, theta_2) -1/(2*pi) * ((l_1(theta_1)-l_1(theta_2)).*l_2_prime(theta_2)-(l_2(theta_1)-l_2(theta_2)).*l_1_prime(theta_2))./(sqrt(l_1_prime(theta_2).^2+l_2_prime(theta_2).^2).*((l_1(theta_1)-l_1(theta_2)).^2+(l_2(theta_1)-l_2(theta_2)).^2));
K_boundary_same_point = @(theta) -1/(4*pi) * (l_2_prime(theta).*l_1_dprime(theta)-l_1_prime(theta).*l_2_dprime(theta))./(l_1_prime(theta).^2+l_2_prime(theta).^2).^(3/2);

K_general = @(x, theta) -1/(2*pi) * ((x(1)-l_1(theta)).*l_2_prime(theta)-(x(2)-l_2(theta)).*l_1_prime(theta))./(sqrt(l_1_prime(theta).^2+l_2_prime(theta).^2).*((x(1)-l_1(theta)).^2+(x(2)-l_2(theta)).^2));

n = 160;
dtheta = 1/n;
theta_mesh = linspace(0, 1-dtheta, n)';
msk_diagonals = logical(diag(ones(n, 1)));

[theta_2, theta_1] = meshgrid(theta_mesh);

A = K_boundary(theta_1, theta_2).*sqrt(l_1_prime(theta_mesh').^2+l_2_prime(theta_mesh').^2)*dtheta;
A(msk_diagonals) = K_boundary_same_point(theta_mesh).*sqrt(l_1_prime(theta_mesh).^2+l_2_prime(theta_mesh).^2)*dtheta + 1/2;

b = f(l_1(theta_mesh), l_2(theta_mesh));
phi_j = A\b;

h = 0.05;
delta = 0.1;
refined_theta_mesh = linspace(0, 1, n*100)';
boundary_X = l_1(refined_theta_mesh);
boundary_Y = l_2(refined_theta_mesh);

R = R_cartesian_mesh_obj(-1.5-h*rand, 1, -1-h*rand, 1, h, boundary_X, boundary_Y);
near_boundary_msk = compute_near_boundary_points(delta, R.R_X, R.R_Y, boundary_X, boundary_Y);

interior_point_idxs = R.R_idxs(R.in_interior & ~near_boundary_msk);
near_boundary_point_idxs = R.R_idxs(R.in_interior & near_boundary_msk);

figure;
scatter(R.R_X(R.in_interior), R.R_Y(R.in_interior))
hold on;
plot(boundary_X(:), boundary_Y(:))

% numeric computation in interior
f_numeric = zeros(size(R.R_X));
for idx = interior_point_idxs'
    f_numeric(idx) = transpose(K_general([R.R_X(idx); R.R_Y(idx)], theta_mesh).*sqrt(l_1_prime(theta_mesh).^2+l_2_prime(theta_mesh).^2).*phi_j) * ones(n, 1) * dtheta;
end

% error calculation for interior point (can do convergence analysis later)
err_interior = max(abs(f(R.R_X(interior_point_idxs), R.R_Y(interior_point_idxs))-f_numeric(interior_point_idxs)))

% near boundary calculation
fc_coeff_phi = fftshift(fft(phi_j.*sqrt(l_1_prime(theta_mesh).^2+l_2_prime(theta_mesh).^2)))/n;


if mod(n, 2) == 0
    freq_mesh = -n/2:(n/2-1);
else
    freq_mesh = (-(n-1)/2):((n-1)/2);
end


exp_term = @(theta, freq) exp(2i*pi*freq*theta);

for idx = near_boundary_point_idxs'
    gk_integral_vals = zeros(n, 1);
    for freq_idx = 1:length(freq_mesh)
        freq = freq_mesh(freq_idx);
        gk_integral_vals(freq_idx) = quadgk(@(theta) K_general([R.R_X(idx); R.R_Y(idx)], theta).*exp_term(theta, freq), 0, 1, 'AbsTol', err_interior);
    end
    f_numeric(idx) = real(fc_coeff_phi' * gk_integral_vals);
end

f_numeric(~R.in_interior) = nan;

figure;
surf(R.R_X, R.R_Y, f_numeric)

final_err = max(abs(f(R.R_X(R.in_interior), R.R_Y(R.in_interior))-f_numeric(R.in_interior)))

function [near_boundary] = compute_near_boundary_points(delta, R_X, R_Y, boundary_X, boundary_Y) 
   %assumes uniform h
   h = R_X(1, 2) - R_X(1, 1);
   x_start = R_X(1, 1);
   y_start = R_Y(1, 1);
   n_x = size(R_X, 2);
   n_y = size(R_X, 1);
   
   near_boundary = false(size(R_X));
   
   straight_threshold = ceil(delta/h)+1;
   diagonal_threshold = ceil(delta/(sqrt(2)*h))+1;
   
   % finds closest cartesian point to each boundary point
   boundary_X_j = maintain_bounds(round((boundary_X-x_start)/h), 0, n_x-1);
   boundary_Y_j = maintain_bounds(round((boundary_Y-y_start)/h), 0, n_y-1);
   
   straight_left_right_idxs = maintain_bounds(boundary_X_j + (-straight_threshold:straight_threshold), 0, n_x-1);
   fixed_Y = boundary_Y_j + zeros(1, straight_threshold*2+1);
   near_boundary(sub2ind(size(near_boundary), fixed_Y(:)+1, straight_left_right_idxs(:)+1)) = true;
   
   straight_up_down_idxs = maintain_bounds(boundary_Y_j + (-straight_threshold:straight_threshold), 0, n_y-1);
   fixed_X = boundary_X_j + zeros(1, straight_threshold*2+1);
   near_boundary(sub2ind(size(near_boundary), straight_up_down_idxs(:)+1, fixed_X(:)+1)) = true;
   
   diag_1_X = maintain_bounds(boundary_X_j + (-diagonal_threshold:diagonal_threshold), 0, n_x-1);
   diag_1_Y = maintain_bounds(boundary_Y_j + (-diagonal_threshold:diagonal_threshold), 0, n_y-1);
   near_boundary(sub2ind(size(near_boundary), diag_1_Y(:)+1, diag_1_X(:)+1)) = true;
   
   diag_2_Y = maintain_bounds(boundary_Y_j + (diagonal_threshold:-1:-diagonal_threshold), 0, n_y-1);
   near_boundary(sub2ind(size(near_boundary), diag_2_Y(:)+1, diag_1_X(:)+1)) = true;
end

function idx_array = maintain_bounds(idx_array, min_bound, max_bound)
   idx_array(idx_array < min_bound) = min_bound;
   idx_array(idx_array > max_bound) = max_bound;
end