clc; clear; close all;

load('data.mat')
u_G = @(x, y) x.^2-y.^2;
p = 4;
G_cf = 1;
M = 10;
int_eps = fc_err/10;
eps_xi_eta = 1e-14;
eps_xy = 1e-14;
n_r = 6;

n_POU = 50;

curve_seq_coarse = Curve_seq_obj();
curr = curve_seq.first_curve;
for i = 1:curve_seq.n_curves
    curve_seq_coarse.add_curve(curr.l_1, curr.l_2, curr.l_1_prime, curr.l_2_prime, curr.l_1_dprime, curr.l_2_dprime, ceil((curr.n-1)*G_cf)+1, nan, nan, nan, nan, nan);
    curr = curr.next_curve;
end

%% Graded Mesh
v = @(s) (1/p-1/2)*(1-2*s).^3+1/p*(2*s-1)+1/2;
v_prime = @(s) 2/p-6*(1/p-1/2)*(1-2*s).^2;

w = @(s) v(s).^p./(v(s).^p+v(1-s).^p);
w_prime = @(s) p*((v_prime(s).*v(s).^(p-1))./(v(s).^p+v(1-s).^p)-(v(s).^(p-1).*v_prime(s)-v(1-s).^(p-1).*v_prime(1-s)).*v(s).^p./(v(s).^p+v(1-s).^p).^2);

%% Interior Patches

[interior_patches, corner_patches_0, corner_patches_1, C2_corner, corner_theta_thresholds,  corner_theta_j_thresholds] =  construct_interior_patches(curve_seq, R.h , M, eps_xi_eta, eps_xy);
figure;
hold on;
for i = 1:curve_seq.n_curves
    [X, Y] = interior_patches{i}.xy_mesh;
    scatter(X(:), Y(:));
    [X, Y] = corner_patches_0{i}.xy_mesh;
    scatter(X(:), Y(:), 'x')
    [X, Y] = corner_patches_1{i}.xy_mesh;
    scatter(X(:), Y(:), 'x')
end
plot(R.boundary_X, R.boundary_Y);

% computing Cartesian points in each patch
in_S_patch_global = false(size(R.R_X));
in_interior_patch = cell(curve_seq.n_curves, 1);
in_corner_patch_0 = cell(curve_seq.n_curves, 1);
in_corner_patch_1 = cell(curve_seq.n_curves, 1);

for i = 1:curve_seq.n_curves
    %S_patch
    interior_patch = interior_patches{i};
    [bound_X, bound_Y] = interior_patch.boundary_mesh_xy(n_r, false);
    in_patch = inpolygon_mesh(R.R_X, R.R_Y, bound_X, bound_Y) & R.in_interior;
    in_interior_patch{i} = in_patch;
    in_S_patch_global = in_S_patch_global | in_patch;
    
    corner_patch_0 = corner_patches_0{i};
    [bound_X, bound_Y] = corner_patch_0.boundary_mesh_xy(n_r, false);
    in_patch = inpolygon_mesh(R.R_X, R.R_Y, bound_X, bound_Y) & R.in_interior;
    in_corner_patch_0{i} = in_patch;
    in_S_patch_global = in_S_patch_global | in_patch;
    
    corner_patch_1 = corner_patches_1{i};
    [bound_X, bound_Y] = corner_patch_1.boundary_mesh_xy(n_r, false);
    in_patch = inpolygon_mesh(R.R_X, R.R_Y, bound_X, bound_Y) & R.in_interior;
    in_corner_patch_1{i} = in_patch;
    in_S_patch_global = in_S_patch_global | in_patch;
end

% Coarse density
R_coarse = 1;

[A_coarse, b_coarse, n_total, curve_n, start_idx, end_idx] = construct_A_b(R_coarse, u_G, curve_seq_coarse, w, w_prime);
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

% disp("started Gauss-Konrod");
% 
% % Gauss-Konrod quadrature for points only in C type patches
% C_nS_idxs = R.R_idxs(in_C_patch_global & ~ in_S_patch_global);
% for idx = C_nS_idxs'
%     u_num_mat(idx) = u_num_near_boundary_global(R.R_X(idx), R.R_Y(idx), gr_phi_coarse, curve_seq_coarse, start_idx, end_idx, curve_n, w, int_eps, M);
% end
% disp("finished Gauss-Konrod")
%%
% evaluating u_num_coarse on cartesian mesh
interior_msk = ~in_S_patch_global & R.in_interior;
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

    sum(to_update, 'all')
    R_fine
end

%%
% evaulating u_num_coarse on interior patch mesh (S patch only)

% refine one more time
R_coarse = R_fine;
R_fine = R_coarse + 1;
u_num_coarse = u_num_fine;

padded_fft_coeffs = [zeros(ceil((n_total*R_fine-n_total)/2), 1); gr_phi_coarse_fft; zeros(floor((n_total*R_fine-n_total)/2), 1)];
gr_phi_fine = R_fine*n_total*real(ifft(ifftshift(padded_fft_coeffs)));

[~, curve_n_fine, start_idx_fine, end_idx_fine] = compute_curve_param(R_fine, curve_seq_coarse);
u_num_fine = @(x, y) u_num_global(x, y, gr_phi_fine, curve_seq_coarse, start_idx_fine, end_idx_fine, curve_n_fine, w);

POU = construct_POU(n_POU);

% excludes boundary
for eta_idx = M:-1:2
    for i = 1:curve_seq.n_curves
        interior_patch = interior_patches{i};
        [patch_X, patch_Y] = interior_patch.xy_mesh;
        for xi_idx = 1:interior_patch.n_xi
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

                [~, curve_n_fine, start_idx_fine, end_idx_fine] = compute_curve_param(R_fine, curve_seq_coarse);
                u_num_fine = @(x, y) u_num_global(x, y, gr_phi_fine, curve_seq_coarse, start_idx_fine, end_idx_fine, curve_n_fine, w);
            end
            interior_patch.f_XY(eta_idx, xi_idx) = u_num_coarse_val;
        end
        
        corner_patch_0 = corner_patches_0{i};
        [patch_X, patch_Y] = corner_patch_0.xy_mesh;
        
        gr_phi_idx_start = start_idx_fine(i);
        gr_phi_idx_end = start_idx_fine(i)+(corner_theta_j_thresholds(i, 1)*R_fine-1)+n_POU-1;
        gr_phi_POU = gr_phi_fine; gr_phi_POU((gr_phi_idx_end-(n_POU-1)+1):gr_phi_idx_end) = (1-POU(1:n_POU-1)).*gr_phi_POU((gr_phi_idx_end-(n_POU-1)+1):gr_phi_idx_end);
        
        for xi_idx = 1:corner_patch_0.n_xi
            corner_patch_0.f_XY(eta_idx, xi_idx) = int_num_segment(patch_X(eta_idx, xi_idx), patch_Y(eta_idx, xi_idx), gr_phi_POU, curve_seq, gr_phi_idx_start, gr_phi_idx_end, start_idx_fine, end_idx_fine, curve_n_fine, w );
        end
        
        corner_patch_1 = corner_patches_1{i};
        [patch_X, patch_Y] = corner_patch_1.xy_mesh;
        gr_phi_idx_start = start_idx_fine(i) + (corner_theta_j_thresholds(i, 2)-1)*R_fine -n_POU+1;
        gr_phi_idx_end = end_idx_fine(i);
        gr_phi_POU = gr_phi_fine; gr_phi_POU(gr_phi_idx_start:gr_phi_idx_start+n_POU-1) = POU.*gr_phi_POU(gr_phi_idx_start:gr_phi_idx_start+n_POU-1);
        
        for xi_idx = 1:corner_patch_1.n_xi
            corner_patch_1.f_XY(eta_idx, xi_idx) = int_num_segment(patch_X(eta_idx, xi_idx), patch_Y(eta_idx, xi_idx), gr_phi_POU, curve_seq, gr_phi_idx_start, gr_phi_idx_end, start_idx_fine, end_idx_fine, curve_n_fine, w );
        end
    end
end

%test
i1 = 1;
i2 = 2;
x = 0.95;
y = 0.03;

figure;

i = i1;
gr_phi_idx_start = start_idx_fine(i) + (corner_theta_j_thresholds(i, 2)-1)*R_fine -n_POU+1;
gr_phi_idx_end = end_idx_fine(i);
POU_interval_1 = gr_phi_idx_start:gr_phi_idx_start+n_POU-1;
end_1 = start_idx_fine(i) + (corner_theta_j_thresholds(i, 2)-1)*R_fine-1;
gr_phi_POU = gr_phi_fine; gr_phi_POU(gr_phi_idx_start:gr_phi_idx_start+n_POU-1) = POU.*gr_phi_POU(gr_phi_idx_start:gr_phi_idx_start+n_POU-1);
i1_contr = int_num_segment(x, y, gr_phi_POU, curve_seq, gr_phi_idx_start, gr_phi_idx_end, start_idx_fine, end_idx_fine, curve_n_fine, w )

mesh = 1:length(gr_phi_POU);

i = i2;
gr_phi_idx_start = start_idx_fine(i);
gr_phi_idx_end = start_idx_fine(i)+(corner_theta_j_thresholds(i, 1)*R_fine-1)+n_POU-1;
POU_interval_0 = (gr_phi_idx_end-(n_POU-1)+1):gr_phi_idx_end;
start_0 = start_idx_fine(i)+(corner_theta_j_thresholds(i, 1)*R_fine);
gr_phi_POU = gr_phi_fine; gr_phi_POU((gr_phi_idx_end-(n_POU-1)+1):gr_phi_idx_end) = (1-POU(1:n_POU-1)).*gr_phi_POU((gr_phi_idx_end-(n_POU-1)+1):gr_phi_idx_end);
i2_contr = int_num_segment(x, y, gr_phi_POU, curve_seq, gr_phi_idx_start, gr_phi_idx_end, start_idx_fine, end_idx_fine, curve_n_fine, w )

gr_phi_POU = gr_phi_fine; 
gr_phi_POU(POU_interval_1) = (1-POU).*gr_phi_POU(POU_interval_1);
gr_phi_POU(POU_interval_0) = POU(1:n_POU-1).*gr_phi_POU(POU_interval_0);
rest_contr = int_num_segment(x, y, gr_phi_POU, curve_seq, start_0, end_1, start_idx_fine, end_idx_fine, curve_n_fine, w)

i1_contr + i2_contr+rest_contr

% boundary value computation
for i = 1:curve_seq.n_curves
    interior_patch = interior_patches{i};
    [patch_X, patch_Y] = interior_patch.xy_mesh;

    interior_patch.f_XY(1, :) = u_G(patch_X(1, corner_theta_j_thresholds(i, 1):corner_theta_j_thresholds(i, 2)), patch_Y(1, corner_theta_j_thresholds(i, 1):corner_theta_j_thresholds(i, 2)));
end

for i = 1:curve_seq.n_curves
    interior_patch = interior_patches{i};

    in_patch = in_interior_patch{i};
    R_patch_idxs = R.R_idxs(in_patch);

    [P_xi, P_eta] = R_xi_eta_inversion(R, interior_patch, in_patch);
    smooth_patch_msk = P_xi <= corner_theta_thresholds(i, 2) & P_xi >= corner_theta_thresholds(i, 1);

    for idx = 1:length(R_patch_idxs)
        if smooth_patch_msk(idx)
            u_num_mat(R_patch_idxs(idx)) = Q_patch_bounded_locally_compute(interior_patch, P_xi(idx), P_eta(idx), M, corner_theta_j_thresholds(i, :));
        end
    end
end
u_num_mat(~R.in_interior) = nan;

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

function int_num = int_num_segment(x, y, gr_phi, curve_seq, gr_phi_idx_start, gr_phi_idx_end, start_idx, end_idx, curve_n, w)
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
    total_n = sum(curve_n);
    if gr_phi_idx_start > gr_phi_idx_end % wrapped around
        total_nodes = gr_phi_idx_end + total_n - gr_phi_idx_start + 1;
    else
        total_nodes = gr_phi_idx_end - gr_phi_idx_start + 1;
    end
    
    curr = curve_seq.first_curve;
    integrating = false;
    int_num = 0;
    curve_i = 1;
    node_i = 0;
    while true
        ds = 1/curve_n(curve_i);
        if integrating  
            curve_int_start_idx = 1;
        elseif ~integrating && gr_phi_idx_start < end_idx(curve_i)
            integrating = true;
            curve_int_start_idx = gr_phi_idx_start - start_idx(curve_i) + 1;
        end
        
        if integrating && total_nodes <= node_i+(curve_n(curve_i)-curve_int_start_idx+1)
            s_mesh = linspace(0, 1-ds, curve_n(curve_i))'; s_mesh = s_mesh(curve_int_start_idx:((gr_phi_idx_end) - start_idx(curve_i)+1));
            theta_mesh = w(s_mesh);

             gr_phi_idxs = (start_idx(curve_i)-1 + (curve_int_start_idx:(gr_phi_idx_end - start_idx(curve_i)+1)));
             scatter3(curr.l_1(theta_mesh), curr.l_2(theta_mesh), gr_phi(gr_phi_idxs));
             hold on;
             int_num = int_num+transpose(K_general([x; y], theta_mesh, curr).*gr_phi(gr_phi_idxs).*sqrt(curr.l_1_prime(theta_mesh).^2+curr.l_2_prime(theta_mesh).^2)) * ones(length(s_mesh), 1) * ds;
             break;
        elseif integrating
            s_mesh = linspace(0, 1-ds, curve_n(curve_i))'; s_mesh = s_mesh(curve_int_start_idx:end);
            theta_mesh = w(s_mesh);

            gr_phi_idxs = (start_idx(curve_i)-1 + (curve_int_start_idx:curve_n(curve_i)));
            scatter3(curr.l_1(theta_mesh), curr.l_2(theta_mesh), gr_phi(gr_phi_idxs));
            hold on;
            int_num = int_num+transpose(K_general([x; y], theta_mesh, curr).*gr_phi(gr_phi_idxs).*sqrt(curr.l_1_prime(theta_mesh).^2+curr.l_2_prime(theta_mesh).^2)) * ones(length(s_mesh), 1) * ds;
            node_i = node_i + curve_n(curve_i);
        end
        
        curr = curr.next_curve;
        curve_i = mod(curve_i, curve_seq.n_curves)+1;
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

function [interior_patches, corner_patches_0, corner_patches_1, C2_corner, corner_theta_thresholds,  corner_theta_j_thresholds] = construct_interior_patches(curve_seq, h_norm, M, eps_xi_eta, eps_xy)
    C2_corner = false(curve_seq.n_curves, 1);
    corner_theta_thresholds = zeros(curve_seq.n_curves, 2);
    corner_theta_j_thresholds =  zeros(curve_seq.n_curves, 2);
    
    curr = curve_seq.first_curve;
    for i = 1:curve_seq.n_curves
        curr_v = [curr.l_1(1); curr.l_2(1)] - [curr.l_1(1-1/(curr.n-1)); curr.l_2(1-1/(curr.n-1))];
        next_v = [curr.next_curve.l_1(1/(curr.next_curve.n-1)); curr.next_curve.l_2(1/(curr.next_curve.n-1))] - [curr.l_1(1); curr.l_2(1)];
        
        % C2-type patch means cross product is positive
        if curr_v(1)*next_v(2) - curr_v(2)*next_v(1) >= 0
            C2_corner(i) = true;
            
            [theta_1, theta_2] = compute_normal_intersection(curr, curr.next_curve, M, h_norm, eps_xy, [1; 0]);
            corner_theta_j_thresholds(i, 2) = floor(theta_1 * (curr.n-1)) + 1;
            corner_theta_thresholds(i, 2) = (corner_theta_j_thresholds(i, 2)-1)/(curr.n-1);
            if i == curve_seq.n_curves
                corner_theta_j_thresholds(1, 1) = ceil(theta_2 * (curve_seq.first_curve.n-1)) + 1;
                corner_theta_thresholds(1, 1) = (corner_theta_j_thresholds(1, 1) -1)/ (curve_seq.first_curve.n-1);
            else
                corner_theta_j_thresholds(i+1, 1) = ceil(theta_2 * (curr.next_curve.n-1)) + 1;
                corner_theta_thresholds(i+1, 1) = (corner_theta_j_thresholds(i+1, 1) -1)/ (curr.next_curve.n-1);
            end
         % C1-type patch otherwise
        else
            corner_theta_j_thresholds(i, 2) = curr.n;
            corner_theta_thresholds(i, 2) = 1;
            if i == curve_seq.n_curves
                corner_theta_j_thresholds(1, 1) =1;
                corner_theta_thresholds(1, 1) = 0;
            else
                corner_theta_j_thresholds(i+1, 1) = 1;
                corner_theta_thresholds(i+1, 1) = 0;
            end
        end
        
        curr = curr.next_curve;
    end
    
    interior_patches = cell(curve_seq.n_curves);
    corner_patches_0 = cell(curve_seq.n_curves);
    corner_patches_1 = cell(curve_seq.n_curves);
    curr = curve_seq.first_curve;
    for i = 1:curve_seq.n_curves
        xi_diff = 1;
        xi_tilde = @(xi) xi;

        nu_norm = @(theta) sqrt(curr.l_1_prime(theta).^2 + curr.l_2_prime(theta).^2);

        M_p_1 = @(xi, eta) curr.l_1(xi_tilde(xi)) - eta.*curr.l_2_prime(xi_tilde(xi))./nu_norm(xi_tilde(xi));
        M_p_2 = @(xi, eta) curr.l_2(xi_tilde(xi)) + eta.*curr.l_1_prime(xi_tilde(xi))./nu_norm(xi_tilde(xi));
        M_p = @(xi, eta) [M_p_1(xi, eta), M_p_2(xi, eta)];

        dM_p_1_dxi = @(xi, eta) xi_diff * (curr.l_1_prime(xi_tilde(xi))-eta*(curr.l_2_dprime(xi_tilde(xi)).*nu_norm(xi_tilde(xi)).^2-curr.l_2_prime(xi_tilde(xi)).*(curr.l_2_dprime(xi_tilde(xi)).*curr.l_2_prime(xi_tilde(xi))+curr.l_1_dprime(xi_tilde(xi)).*curr.l_1_prime(xi_tilde(xi))))./nu_norm(xi_tilde(xi)).^3);
        dM_p_2_dxi = @(xi, eta) xi_diff * (curr.l_2_prime(xi_tilde(xi))+eta*(curr.l_1_dprime(xi_tilde(xi)).*nu_norm(xi_tilde(xi)).^2-curr.l_1_prime(xi_tilde(xi)).*(curr.l_2_dprime(xi_tilde(xi)).*curr.l_2_prime(xi_tilde(xi))+curr.l_1_dprime(xi_tilde(xi)).*curr.l_1_prime(xi_tilde(xi))))./nu_norm(xi_tilde(xi)).^3);
        dM_p_1_deta = @(xi, eta) -curr.l_2_prime(xi_tilde(xi)) ./ nu_norm(xi_tilde(xi));
        dM_p_2_deta = @(xi, eta) curr.l_1_prime(xi_tilde(xi)) ./ nu_norm(xi_tilde(xi));

        J = @(v) [dM_p_1_dxi(v(1), v(2)), dM_p_1_deta(v(1), v(2)); dM_p_2_dxi(v(1), v(2)), dM_p_2_deta(v(1), v(2))];
        
        n_xi_interior = corner_theta_j_thresholds(i, 2) - corner_theta_j_thresholds(i, 1)+1;
        interior_patches{i} = Q_patch_obj(M_p, J, eps_xi_eta, eps_xy, n_xi_interior, M, corner_theta_thresholds(i, 1), corner_theta_thresholds(i, 2), 0, (M-1)*h_norm, zeros(M, n_xi_interior));

        if corner_theta_thresholds(i, 1) ~= 0
            n_xi_corner =  corner_theta_j_thresholds(i, 1);
            corner_patches_0{i} = Q_patch_obj(M_p, J, eps_xi_eta, eps_xy, n_xi_corner, M, 0, corner_theta_thresholds(i, 1), 0, (M-1)*h_norm, zeros(M, n_xi_corner));
        else
            corner_patches_0{i} = nan;
        end
        if corner_theta_thresholds(i, 2) ~= 1
            n_xi_corner = curr.n - corner_theta_j_thresholds(i, 2) + 1;
            corner_patches_1{i} = Q_patch_obj(M_p, J, eps_xi_eta, eps_xy, n_xi_corner, M, corner_theta_thresholds(i, 2), 1, 0, (M-1)*h_norm, zeros(M, n_xi_corner));
        else
            corner_patches_1{i} = nan;
        end
        
        curr = curr.next_curve;
    end
end

function  [theta_1, theta_2] = compute_normal_intersection(curve_1, curve_2, M, h_norm, eps_xy, initial_guess)
    nu_norm = @(theta, curve) sqrt(curve.l_1_prime(theta).^2 + curve.l_2_prime(theta).^2);
    err_guess_x = @(theta_1, theta_2) curve_1.l_1(theta_1) - (M-1)*h_norm*(curve_1.l_2_prime(theta_1)./nu_norm(theta_1, curve_1)) - (curve_2.l_1(theta_2) - (M-1)*h_norm*(curve_2.l_2_prime(theta_2)./nu_norm(theta_2, curve_2)));
    err_guess_y = @(theta_1, theta_2) curve_1.l_2(theta_1) + (M-1)*h_norm*(curve_1.l_1_prime(theta_1)./nu_norm(theta_1, curve_1)) - (curve_2.l_2(theta_2) + (M-1)*h_norm*(curve_2.l_1_prime(theta_2)./nu_norm(theta_2, curve_2)));

    derr_x_d1 = @(theta_1, theta_2) curve_1.l_1_prime(theta_1)-(M-1)*h_norm*(curve_1.l_2_dprime(theta_1).*nu_norm(theta_1, curve_1).^2-curve_1.l_2_prime(theta_1).*(curve_1.l_2_dprime(theta_1).*curve_1.l_2_prime(theta_1)+curve_1.l_1_dprime(theta_1).*curve_1.l_1_prime(theta_1)))./nu_norm(theta_1, curve_1).^3;
    derr_y_d1 = @(theta_1, theta_2) curve_1.l_2_prime(theta_1)+(M-1)*h_norm*(curve_1.l_1_dprime(theta_1).*nu_norm(theta_1, curve_1).^2-curve_1.l_1_prime(theta_1).*(curve_1.l_2_dprime(theta_1).*curve_1.l_2_prime(theta_1)+curve_1.l_1_dprime(theta_1).*curve_1.l_1_prime(theta_1)))./nu_norm(theta_1, curve_1).^3;
    derr_x_d2 = @(theta_1, theta_2) -(curve_2.l_1_prime(theta_2)-(M-1)*h_norm*(curve_2.l_2_dprime(theta_2).*nu_norm(theta_2, curve_2).^2-curve_2.l_2_prime(theta_2).*(curve_2.l_2_dprime(theta_2).*curve_2.l_2_prime(theta_2)+curve_2.l_1_dprime(theta_2).*curve_2.l_1_prime(theta_2)))./nu_norm(theta_2, curve_2).^3);
    derr_y_d2 = @(theta_1, theta_2) -( curve_2.l_2_prime(theta_2)+(M-1)*h_norm*(curve_2.l_1_dprime(theta_2).*nu_norm(theta_2, curve_2).^2-curve_2.l_1_prime(theta_2).*(curve_2.l_2_dprime(theta_2).*curve_2.l_2_prime(theta_2)+curve_2.l_1_dprime(theta_2).*curve_2.l_1_prime(theta_2)))./nu_norm(theta_2, curve_2).^3);

    err_guess = @(v) [err_guess_x(v(1), v(2)); err_guess_y(v(1), v(2))];
    J_err = @(v) [derr_x_d1(v(1), v(2)) derr_x_d2(v(1), v(2)); derr_y_d1(v(1), v(2)) derr_y_d2(v(1), v(2))];

    [v_guess, converged] = newton_solve(err_guess, J_err, initial_guess, eps_xy, 100);
    theta_1 = v_guess(1); theta_2 = v_guess(2);
    if ~converged
        warning('Nonconvergence in normal-boundary intersection')
    end
end

function [f_xy, in_range] = Q_patch_bounded_locally_compute(Q_patch, xi, eta, M, corner_xi_j_thresholds)
    % locally_compute uses two step one dimensional polynomial interpolation to estimate the
    %   value of f for any xi, eta within the bounds of the patch (not necessarily on discrete mesh)
    %
    % Input parameters:
    %    xi (double) - scalar xi value within bounds of patch in
    %       parameter space
    %    eta (double) - scalar eta value within bounds of patch in
    %       parameter space
    %   M (int) - number of interpolation points used for each one
    %       dimensional procedure

    % check if xi, eta are within bounds of Q
    if ~Q_patch.in_patch(xi, eta)
        f_xy = nan;
        in_range = false;
        warning('locally computing for point not in patch');
        return
    end

    in_range = true;
    [h_xi, h_eta] = Q_patch.h_mesh();

    % j to the immediate left of point
    xi_j = floor((xi-Q_patch.xi_start)/h_xi);
    eta_j = floor((eta-Q_patch.eta_start)/h_eta);

    half_M =  floor(M/2);
    if mod(M, 2) ~= 0
        interpol_xi_j_mesh = transpose(xi_j-half_M:xi_j+half_M);
        interpol_eta_j_mesh = transpose(eta_j-half_M:eta_j+half_M);
    else
        interpol_xi_j_mesh = transpose(xi_j-half_M+1:xi_j+half_M);
        interpol_eta_j_mesh = transpose(eta_j-half_M+1:eta_j+half_M);
    end

    interpol_xi_j_mesh = shift_idx_mesh(interpol_xi_j_mesh, corner_xi_j_thresholds(1)-1, corner_xi_j_thresholds(2)-1);
    interpol_eta_j_mesh = shift_idx_mesh(interpol_eta_j_mesh, 0, Q_patch.n_eta-1);

    interpol_xi_mesh = h_xi*interpol_xi_j_mesh + Q_patch.xi_start;
    interpol_eta_mesh = h_eta*interpol_eta_j_mesh + Q_patch.eta_start;

    % first 1D interpolation
    interpol_xi_exact = zeros(M, 1);
    for horz_idx = 1:M
        mu = [mean(interpol_xi_mesh), std(interpol_xi_mesh)];
        interpol_val = Q_patch.f_XY(interpol_eta_j_mesh(horz_idx)+1, interpol_xi_j_mesh+1)';
        interpol_xi_exact(horz_idx) = barylag([(interpol_xi_mesh-mu(1))/mu(2), interpol_val], (xi-mu(1))/mu(2));
    end
     % second 1D interpolation
    f_xy = barylag([interpol_eta_mesh, interpol_xi_exact], eta);
end

function POU = construct_POU(n_POU)
    w_1D = @(x) erfc(6*(-2*x+1))/2;
    w_1 = @(s) w_1D(s);
    w_2 = @(s) w_1D(-s+1);
    POU_mesh = linspace(0, 1, n_POU)';
    
    POU = w_1(POU_mesh)./(w_1(POU_mesh)+w_2(POU_mesh));
end