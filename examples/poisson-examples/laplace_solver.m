clc; clear; close all;
load('geo_data_teardrop.mat')

coarse_factor = 1/10;
curve_seq_coarse = Curve_seq_obj();
curr = curve_seq.first_curve;
for i = 1:curve_seq.n_curves
    curve_seq_coarse.add_curve(curr.l_1, curr.l_2, curr.l_1_prime, curr.l_2_prime, curr.l_1_dprime, curr.l_2_dprime, ceil((curr.n-1)*coarse_factor)+1, nan, nan, nan, nan, nan);
    curr = curr.next_curve;
end

%% Graded Mesh
p = 2;
f = @(x, y) x.^2 - y.^2;

v = @(s) (1/p-1/2)*(1-2*s).^3+1/p*(2*s-1)+1/2;
v_prime = @(s) 2/p-6*(1/p-1/2)*(1-2*s).^2;

w = @(s) v(s).^p./(v(s).^p+v(1-s).^p);
w_prime = @(s) p*((v_prime(s).*v(s).^(p-1))./(v(s).^p+v(1-s).^p)-(v(s).^(p-1).*v_prime(s)-v(1-s).^(p-1).*v_prime(1-s)).*v(s).^p./(v(s).^p+v(1-s).^p).^2);

%% Interior Patches
M = 7;
curr = curve_seq.first_curve;
for i=1:curve_seq.n_curves
    curr.h_norm = R.h;
    curr = curr.next_curve;
end
interior_patches = curve_seq.construct_patches(@(x, y) ones(size(x)), M, 1e-13, 1e-13);

% computing Cartesian points in each patch
in_S_patch_global = false(size(R.R_X));
in_C_patch_global = false(size(R.R_X));

in_S_patch = cell(curve_seq.n_curves, 1);
for i = 1:curve_seq.n_curves
    %S_patch
    S_interior_patch = interior_patches{2*i-1}.Q;
    [bound_X, bound_Y] = S_interior_patch.boundary_mesh_xy(false);
    in_patch = inpolygon_mesh(R.R_X, R.R_Y, bound_X, bound_Y) & R.in_interior;
    in_S_patch{i} = in_patch;
    in_S_patch_global = in_S_patch_global | in_patch;
    
    %C_patch 
    C_L_interior_patch = interior_patches{2*i}.L;
    [bound_X, bound_Y] = C_L_interior_patch.boundary_mesh_xy(false);
    in_patch = inpolygon_mesh(R.R_X, R.R_Y, bound_X, bound_Y) & R.in_interior;
    in_C_patch_global = in_C_patch_global | in_patch;
    
    C_W_interior_patch = interior_patches{2*i}.W;
    [bound_X, bound_Y] = C_W_interior_patch.boundary_mesh_xy(false);
    in_patch = inpolygon_mesh(R.R_X, R.R_Y, bound_X, bound_Y) & R.in_interior;
    in_C_patch_global = in_C_patch_global | in_patch;
end

%% 
% Integration epsilon
int_eps = 1e-6;

% Coarse density
R_coarse = 1;

[A_coarse, b_coarse, n_total, curve_n, start_idx, end_idx] = construct_A_b(R_coarse, f, curve_seq_coarse, w, w_prime);
phi_coarse = A_coarse\b_coarse;

gr_phi_coarse = zeros(n_total, 1);
curr = curve_seq_coarse.first_curve;
for i = 1:curve_seq.n_curves
    ds = 1/curve_n(i);
    s_mesh = linspace(0, 1-ds, curve_n(i))';

    gr_phi_coarse(start_idx(i):end_idx(i)) = w_prime(s_mesh).*phi_coarse(start_idx(i):end_idx(i));
    curr = curr.next_curve;
end
gr_phi_coarse_fft = fftshift(fft(gr_phi_coarse))/n_total;

u_num_mat = zeros(R.n_y, R.n_x);

% coarse density u_num
u_num_coarse = @(x, y) u_num_global(x, y, gr_phi_coarse, curve_seq_coarse, start_idx, end_idx, curve_n, w);

% Gauss-Konrod quadrature for points only in C type patches
C_nS_idxs = R.R_idxs(in_C_patch_global & ~ in_S_patch_global);
for idx = C_nS_idxs'
    u_num_mat(idx) = u_num_near_boundary_global(R.R_X(idx), R.R_Y(idx), gr_phi_coarse_fft, curve_seq_coarse, w, int_eps);
end
%%
% evaluating u_num_coarse on cartesian mesh
interior_msk = ~in_C_patch_global & ~in_S_patch_global & R.in_interior;
interior_idxs = R.R_idxs(interior_msk);
for idx = interior_idxs'
    u_num_mat(idx) = u_num_coarse(R.R_X(idx), R.R_Y(idx));
end

%%
% finer density via fft and ifft
to_update = interior_msk;
u_num_mat_fine = zeros(R.n_y, R.n_x);
R_fine = 2;

while sum(to_update, 'all') > 0
    padded_fft_coeffs = [zeros(ceil((n_total*R_fine-n_total)/2), 1); gr_phi_coarse_fft; zeros(floor((n_total*R_fine-n_total)/2), 1)];
    gr_phi_fine = R_fine*n_total*real(ifft(ifftshift(padded_fft_coeffs)));

    [~, curve_n_fine, start_idx_fine, end_idx_fine] = compute_curve_param(R_fine, curve_seq_coarse);
    u_num_fine = @(x, y) u_num_global(x, y, gr_phi_fine, curve_seq_coarse, start_idx_fine, end_idx_fine, curve_n_fine, w);

    for idx = R.R_idxs(to_update)'
        u_num_mat_fine(idx) = u_num_fine(R.R_X(idx), R.R_Y(idx));
    end

    to_update = to_update & abs(u_num_mat_fine - u_num_mat) > int_eps;
    u_num_mat(to_update) = u_num_mat_fine(to_update);
    R_fine = R_fine + 1;
end

%%
% evaulating u_num_coarse on interior patch mesh (S patch only)

% refine one more time
R_coarse = R_fine;
R_fine = R_coarse + 1;
u_num_coarse = u_num_fine;

padded_fft_coeffs = [zeros(ceil((n_total*R_fine-n_total)/2), 1); gr_phi_coarse_fft; zeros(floor((n_total*R_fine-n_total)/2), 1)];
gr_phi_fine = R_fine*n_total*real(ifft(ifftshift(padded_fft_coeffs)));

[n_total_fine, curve_n_fine, start_idx_fine, end_idx_fine] = compute_curve_param(R_fine, curve_seq_coarse);
u_num_fine = @(x, y) u_num_global(x, y, gr_phi_fine, curve_seq_coarse, start_idx_fine, end_idx_fine, curve_n_fine, w);

% excludes boundary
for eta_idx = M:-1:2
    for i = 1:curve_seq.n_curves
        Q_patch = interior_patches{2*i-1}.Q;
        [patch_X, patch_Y] = Q_patch.xy_mesh;
        for xi_idx = 1:size(patch_X, 2)
            while true
                u_num_coarse_val = u_num_coarse(patch_X(eta_idx, xi_idx), patch_Y(eta_idx, xi_idx));
                u_num_fine_val = u_num_fine(patch_X(eta_idx, xi_idx), patch_Y(eta_idx, xi_idx));
                
                if abs(u_num_coarse_val - u_num_fine_val) < int_eps
                    break;
                end
                R_coarse = R_fine;
                R_fine = R_coarse + 1;
                u_num_coarse = u_num_fine;

                padded_fft_coeffs = [zeros(ceil((n_total*R_fine-n_total)/2), 1); gr_phi_coarse_fft; zeros(floor((n_total*R_fine-n_total)/2), 1)];
                gr_phi_fine = R_fine*n_total*real(ifft(ifftshift(padded_fft_coeffs)));

                [n_total_fine, curve_n_fine, start_idx_fine, end_idx_fine] = compute_curve_param(R_fine, curve_seq_coarse);
                u_num_fine = @(x, y) u_num_global(x, y, gr_phi_fine, curve_seq_coarse, start_idx_fine, end_idx_fine, curve_n_fine, w);
            end
            Q_patch.f_XY(eta_idx, xi_idx) = u_num_coarse(patch_X(eta_idx, xi_idx), patch_Y(eta_idx, xi_idx));
        end
    end
end

% boundary value computation
for i = 1:curve_seq.n_curves
    Q_patch = interior_patches{2*i-1}.Q;
    [patch_X, patch_Y] = Q_patch.xy_mesh;
    
    Q_patch.f_XY(1, :) = f(patch_X(1, :), patch_Y(1, :));
end

for i = 1:curve_seq.n_curves
    Q_patch = interior_patches{2*i-1}.Q;
    
    [bound_X, bound_Y] = Q_patch.boundary_mesh_xy(false);
    in_patch = in_S_patch{i};
    R_patch_idxs = R.R_idxs(in_patch);

    [P_xi, P_eta] = R_xi_eta_inversion(R, Q_patch, in_patch);
    
    for idx = 1:length(R_patch_idxs)
        u_num_mat(R_patch_idxs(idx)) = Q_patch.locally_compute(P_xi(idx), P_eta(idx), M);
    end
 end

function [n_total, curve_n, start_idx, end_idx] = compute_curve_param(R, curve_seq)
    curve_n = zeros(curve_seq.n_curves, 1);
    curr = curve_seq.first_curve;
    for i = 1:curve_seq.n_curves
        curve_n(i) = ceil((curr.n-1)*R);
        curr = curr.next_curve;
    end

    n_total = sum(curve_n);
    start_idx = cumsum([1, curve_n(1:end-1)'])';
    end_idx = start_idx+curve_n-1;
end

function [A, b, n_total, curve_n, start_idx, end_idx] = construct_A_b(R, f, curve_seq, w, w_prime)
    [n_total, curve_n, start_idx, end_idx] = compute_curve_param(R, curve_seq);

    K_boundary = @(theta_1, theta_2, curve_1, curve_2) ...
        -1 / (2 * pi) * ( ...
            (curve_1.l_1(theta_1) - curve_2.l_1(theta_2)) .* curve_2.l_2_prime(theta_2) - ...
            (curve_1.l_2(theta_1) - curve_2.l_2(theta_2)) .* curve_2.l_1_prime(theta_2) ...
        ) ./ ( ...
            sqrt(curve_2.l_1_prime(theta_2).^2 + curve_2.l_2_prime(theta_2).^2) .* ...
            ( ...
                (curve_1.l_1(theta_1) - curve_2.l_1(theta_2)).^2 + ...
                (curve_1.l_2(theta_1) - curve_2.l_2(theta_2)).^2 ...
            ) ...
        );

    K_boundary_same_point = @(theta, curve) ...
        -1 / (4 * pi) * ( ...
            curve.l_2_prime(theta) .* curve.l_1_dprime(theta) - ...
            curve.l_1_prime(theta) .* curve.l_2_dprime(theta) ...
        ) ./ ( ...
            (curve.l_1_prime(theta).^2 + curve.l_2_prime(theta).^2).^(3/2) ...
        );

    A = zeros(n_total, n_total);
    b = zeros(n_total, 1);


    curr_y = curve_seq.first_curve;
    for i = 1:curve_seq.n_curves
        ds_y = 1/curve_n(i);
        s_mesh_y = linspace(0, 1-ds_y, curve_n(i))';
        theta_mesh_y = w(s_mesh_y);
        b(start_idx(i):end_idx(i)) = f(curr_y.l_1(theta_mesh_y), curr_y.l_2(theta_mesh_y));

        curr_x = curve_seq.first_curve;
        for j = 1:curve_seq.n_curves
            ds_x = 1/curve_n(j);
            s_mesh_x = linspace(0, 1-ds_x, curve_n(j))';
            theta_mesh_x = w(s_mesh_x);
            [theta_2, theta_1] = meshgrid(theta_mesh_y, theta_mesh_x);

            A_local = w_prime(s_mesh_y').*K_boundary(theta_1, theta_2, curr_x, curr_y).*sqrt(curr_y.l_1_prime(theta_mesh_y').^2+curr_y.l_2_prime(theta_mesh_y').^2)*ds_y;

            if i == j
                msk_diagonals = logical(diag(ones(curve_n(j), 1)));
                A_local(msk_diagonals) = w_prime(s_mesh_y).*K_boundary_same_point(theta_mesh_y, curr_x).*sqrt(curr_x.l_1_prime(theta_mesh_x).^2+curr_x.l_2_prime(theta_mesh_x).^2)*ds_x + 1/2;
            end

            A(start_idx(j):end_idx(j), start_idx(i):end_idx(i)) = A_local;

            curr_x = curr_x.next_curve;
        end

        curr_y = curr_y.next_curve;
    end
end

function u_num_b = u_num_b(curve_idx, s_idx, gr_phi, curve_seq, start_idx, end_idx, curve_n, w, w_prime)
     K_boundary = @(theta_1, theta_2, curve_1, curve_2) ...
        -1 / (2 * pi) * ( ...
            (curve_1.l_1(theta_1) - curve_2.l_1(theta_2)) .* curve_2.l_2_prime(theta_2) - ...
            (curve_1.l_2(theta_1) - curve_2.l_2(theta_2)) .* curve_2.l_1_prime(theta_2) ...
        ) ./ ( ...
            sqrt(curve_2.l_1_prime(theta_2).^2 + curve_2.l_2_prime(theta_2).^2) .* ...
            ( ...
                (curve_1.l_1(theta_1) - curve_2.l_1(theta_2)).^2 + ...
                (curve_1.l_2(theta_1) - curve_2.l_2(theta_2)).^2 ...
            ) ...
        );

    K_boundary_same_point = @(theta, curve) ...
        -1 / (4 * pi) * ( ...
            curve.l_2_prime(theta) .* curve.l_1_dprime(theta) - ...
            curve.l_1_prime(theta) .* curve.l_2_dprime(theta) ...
        ) ./ ( ...
            (curve.l_1_prime(theta).^2 + curve.l_2_prime(theta).^2).^(3/2) ...
        );
    
    compute_curve = nan;
    curr = curve_seq.first_curve;
    for i = 1:size(start_idx)
        if i == curve_idx
            compute_curve = curr;
            compute_s = (s_idx-1)/curve_n(i);
            compute_theta = w(compute_s);
            break;
        end
        curr = curr.next_curve;
    end
    
    curr = curve_seq.first_curve;
    u_num_b = 0;
    for i = 1:size(start_idx)
        ds = 1/curve_n(i);
        s_mesh = linspace(0, 1-ds, curve_n(i))';
        theta_mesh = w(s_mesh);

        int_vec =  K_boundary(compute_theta, theta_mesh, compute_curve, curr).*sqrt(curr.l_1_prime(theta_mesh).^2+curr.l_2_prime(theta_mesh).^2)*ds;
        
        if curve_idx == i
            int_vec(s_idx) = K_boundary_same_point(compute_theta, compute_curve).*sqrt(compute_curve.l_1_prime(compute_theta).^2+compute_curve.l_2_prime(compute_theta).^2)*ds;
            u_num_b = u_num_b + 1/2*gr_phi(start_idx(curve_idx)+s_idx-1)./w_prime(compute_s);
        end

        u_num_b = u_num_b + int_vec' * gr_phi(start_idx(i):end_idx(i));
        curr = curr.next_curve;
    end
end


function u_num = u_num_global(x, y, gr_phi, curve_seq, start_idx, end_idx, curve_n, w)
    K_general = @(x, theta, curve) ...
    -1 / (2 * pi) * ( ...
        (x(1) - curve.l_1(theta)) .* curve.l_2_prime(theta) - ...
        (x(2) - curve.l_2(theta)) .* curve.l_1_prime(theta) ...
    ) ./ ( ...
        sqrt(curve.l_1_prime(theta).^2 + curve.l_2_prime(theta).^2) .* ...
        ( ...
            (x(1) - curve.l_1(theta)).^2 + ...
            (x(2) - curve.l_2(theta)).^2 ...
        ) ...
    );

   curr = curve_seq.first_curve;
   u_num = 0;
   for i = 1:size(start_idx)
        ds = 1/curve_n(i);
        s_mesh = linspace(0, 1-ds, curve_n(i))';
        theta_mesh = w(s_mesh);

        u_num =  u_num + transpose(K_general([x; y], theta_mesh, curr).*gr_phi(start_idx(i):end_idx(i)).*sqrt(curr.l_1_prime(theta_mesh).^2+curr.l_2_prime(theta_mesh).^2)) * ones(curve_n(i), 1) * ds;
        curr = curr.next_curve;
    end
end

function u_num = u_num_near_boundary_global(x, y, gr_phi_fft, curve_seq, w, int_eps)
    K_general = @(x, theta, curve) ...
    -1 / (2 * pi) * ( ...
        (x(1) - curve.l_1(theta)) .* curve.l_2_prime(theta) - ...
        (x(2) - curve.l_2(theta)) .* curve.l_1_prime(theta) ...
    ) ./ ( ...
        sqrt(curve.l_1_prime(theta).^2 + curve.l_2_prime(theta).^2) .* ...
        ( ...
            (x(1) - curve.l_1(theta)).^2 + ...
            (x(2) - curve.l_2(theta)).^2 ...
        ) ...
    );

    n = length(gr_phi_fft);

    if mod(n, 2) == 0
        freq_mesh = -n/2:(n/2-1);
    else
        freq_mesh = (-(n-1)/2):((n-1)/2);
    end

    exp_term = @(s, freq) exp(2i*pi*freq*s);

    gk_integral_vals = zeros(n, 1);
    curr = curve_seq.first_curve;
    for freq_idx = 1:length(freq_mesh)
        freq = freq_mesh(freq_idx);
        
        for i = 1:curve_seq.n_curves
            gk_integral_vals(freq_idx) = gk_integral_vals(freq_idx) + quadgk(@(s) K_general([x; y], w(s), curr).*exp_term(s, freq).*sqrt(curr.l_1_prime(w(s)).^2+curr.l_2_prime(w(s)).^2), 0, 1, 'AbsTol', max(int_eps/n, 1e-14));
            curr = curr.next_curve;
        end
    end
    u_num = real(gr_phi_fft' * gk_integral_vals);
end