clear; close all;

f = @(x, y) x.^2-y.^2;

% teardrop
l_1 = @(theta) 2*sin(theta*pi);
l_2 = @(theta) -sin(theta*2*pi);
l_1_prime = @(theta) 2*pi*cos(theta*pi);
l_2_prime = @(theta) -2*pi*cos(theta*2*pi);
l_1_dprime = @(theta) -2*pi^2*sin(theta*pi);
l_2_dprime = @(theta) 4*pi^2*sin(theta*2*pi);

K_boundary = @(theta_1, theta_2) -1/(2*pi) * ((l_1(theta_1)-l_1(theta_2)).*l_2_prime(theta_2)-(l_2(theta_1)-l_2(theta_2)).*l_1_prime(theta_2))./(sqrt(l_1_prime(theta_2).^2+l_2_prime(theta_2).^2).*((l_1(theta_1)-l_1(theta_2)).^2+(l_2(theta_1)-l_2(theta_2)).^2));
K_boundary_same_point = @(theta) -1/(4*pi) * (l_2_prime(theta).*l_1_dprime(theta)-l_1_prime(theta).*l_2_dprime(theta))./(l_1_prime(theta).^2+l_2_prime(theta).^2).^(3/2);

K_general = @(x, theta) -1/(2*pi) * ((x(1)-l_1(theta)).*l_2_prime(theta)-(x(2)-l_2(theta)).*l_1_prime(theta))./(sqrt(l_1_prime(theta).^2+l_2_prime(theta).^2).*((x(1)-l_1(theta)).^2+(x(2)-l_2(theta)).^2));

p = 20;
v = @(s) (1/p-1/2)*(1-2*s).^3+1/p*(2*s-1)+1/2;
v_prime = @(s) 2/p-6*(1/p-1/2)*(1-2*s).^2;

w = @(s) v(s).^p./(v(s).^p+v(1-s).^p);
w_prime = @(s) p*((v_prime(s).*v(s).^(p-1))./(v(s).^p+v(1-s).^p)-(v(s).^(p-1).*v_prime(s)-v(1-s).^(p-1).*v_prime(1-s)).*v(s).^p./(v(s).^p+v(1-s).^p).^2);

n = 320;
ds = 1/n;
s_mesh = linspace(0, 1-ds, n)';
theta_mesh = w(s_mesh);

figure;
scatter(l_1(theta_mesh), l_2(theta_mesh));

msk_diagonals = logical(diag(ones(n, 1)));

[s_2, s_1] = meshgrid(s_mesh);
[theta_2, theta_1] = meshgrid(theta_mesh);

A = w_prime(s_mesh').*K_boundary(theta_1, theta_2).*sqrt(l_1_prime(theta_mesh').^2+l_2_prime(theta_mesh').^2)*ds;
A(msk_diagonals) = w_prime(s_mesh).*K_boundary_same_point(theta_mesh).*sqrt(l_1_prime(theta_mesh).^2+l_2_prime(theta_mesh).^2)*ds + 1/2;

b = f(l_1(theta_mesh), l_2(theta_mesh));
phi_j = A\b;

fc_coeff_phi = fftshift(fft(w_prime(s_mesh).*phi_j.*sqrt(l_1_prime(theta_mesh).^2+l_2_prime(theta_mesh).^2)))/n;

R_fac = 2;
n_R = n*R_fac;
ds_R = ds/R_fac;
s_R_mesh = linspace(0, 1-ds_R, n_R)';
theta_R_mesh = w(s_R_mesh);

padded_fc_coeffs = [zeros(floor((n_R-n)/2), 1); fc_coeff_phi; zeros(ceil((n_R-n)/2), 1)];
phi_j_R = n_R*real(ifft(ifftshift(padded_fc_coeffs)))./w_prime(s_R_mesh)./sqrt(l_1_prime(theta_R_mesh).^2+l_2_prime(theta_R_mesh).^2); phi_j_R(1) = phi_j(1);

u_num_R = @(x, y)  transpose(w_prime(s_R_mesh).*K_general([x; y], theta_R_mesh).*sqrt(l_1_prime(theta_R_mesh).^2+l_2_prime(theta_R_mesh).^2).*phi_j_R) * ones(n_R, 1) * ds_R;

h = 0.05;
delta = h;
refined_theta_mesh = linspace(0, 1, n*100)';
boundary_X = l_1(refined_theta_mesh);
boundary_Y = l_2(refined_theta_mesh);

R = R_cartesian_mesh_obj(0-h*rand, 2, -1-h*rand, 1, h, boundary_X, boundary_Y);
near_boundary_msk = compute_near_boundary_points(delta, R.R_X, R.R_Y, boundary_X, boundary_Y);

interior_point_idxs = R.R_idxs(R.in_interior & ~near_boundary_msk);
near_boundary_point_idxs = R.R_idxs(R.in_interior & near_boundary_msk);
length(near_boundary_point_idxs)
figure;
scatter(R.R_X(R.in_interior), R.R_Y(R.in_interior))
hold on;
plot(boundary_X(:), boundary_Y(:))

% numeric computation in interior
f_numeric = zeros(size(R.R_X));
for idx = interior_point_idxs'
    f_numeric(idx) = u_num_R(R.R_X(idx), R.R_Y(idx));
end

% error calculation for interior point (can do convergence analysis later)
err_interior = max(abs(f(R.R_X(interior_point_idxs), R.R_Y(interior_point_idxs))-f_numeric(interior_point_idxs)))

figure;
scatter3(R.R_X(interior_point_idxs), R.R_Y(interior_point_idxs), f_numeric(interior_point_idxs))
hold on;
scatter3(R.R_X(interior_point_idxs), R.R_Y(interior_point_idxs), f(R.R_X(interior_point_idxs), R.R_Y(interior_point_idxs)))

% near boundary calculation

if mod(n, 2) == 0
    freq_mesh = -n/2:(n/2-1);
else
    freq_mesh = (-(n-1)/2):((n-1)/2);
end


exp_term = @(s, freq) exp(2i*pi*freq*s);

for idx = near_boundary_point_idxs'
    gk_integral_vals = zeros(n, 1);
    for freq_idx = 1:length(freq_mesh)
        freq = freq_mesh(freq_idx);
        gk_integral_vals(freq_idx) = quadgk(@(s) K_general([R.R_X(idx); R.R_Y(idx)], w(s)).*exp_term(s, freq), 0, 1, 'AbsTol', err_interior);
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