function [u_num_mat] = laplace_solver_new(R, curve_seq, u_G, cfac, p, M, int_eps, eps_xi_eta, eps_xy, rho, rho_intp)
%LAPLACE_SOLVER_NEW Summary of this function goes here
%   Detailed explanation goes here

%% Computes the n discretization values
curve_n_rho1 = zeros(curve_seq.n_curves, 1);
curr = curve_seq.first_curve;
for i = 1:curve_seq.n_curves
    curve_n_rho1(i) = ceil((curr.n-1)*cfac);
    curr = curr.next_curve;
end
curve_param_rho1 = curve_param_obj(curve_n_rho1);

% Constructs U
n_POU = 50;
U = construct_POU(n_POU);

%% Constructs IE curve seq obj

IE_curve_seq = IE_curve_seq_obj(curve_seq, p);
[A_rho1, b_rho1] = IE_curve_seq.construct_A_b(curve_param_rho1, u_G);
gr_phi_rho1 = A_rho1 \ b_rho1;
gr_phi_fft_rho1 = fftshift(fft(gr_phi_rho1))/curve_param_rho1.n_total;

[s_patches, c_0_patches, c_1_patches] = IE_curve_seq.construct_interior_patches(curve_param_rho1, R.h/rho_intp, M, eps_xi_eta, eps_xy);
[well_interior_msk, s_patch_msks, c_0_patch_msks, c_1_patch_msks] = gen_R_msks(R, rho, s_patches, c_0_patches, c_1_patches);

%% Fills in points that are well within the domain
u_num_mat = zeros(size(R.f_R));
u_num_mat_fine = zeros(size(R.f_R));

for idx = R.R_idxs(well_interior_msk)'
    u_num_mat(idx) = IE_curve_seq.u_num(R.R_X(idx), R.R_Y(idx), curve_param_rho1, gr_phi_rho1);
end

to_update = well_interior_msk;
rho_fine_IE = 2;

while sum(to_update, 'all') > 0
    [gr_phi_fine, curve_param_fine] = refine_gr_phi(curve_param_rho1, rho_fine_IE, gr_phi_fft_rho1);

    for idx = R.R_idxs(to_update)'
        u_num_mat_fine(idx) = IE_curve_seq.u_num(R.R_X(idx), R.R_Y(idx), curve_param_fine, gr_phi_fine);
    end

    to_update = to_update & abs(u_num_mat_fine - u_num_mat) > int_eps;
    u_num_mat(to_update) = u_num_mat_fine(to_update);
    rho_fine_IE = rho_fine_IE + 1;

    sum(to_update, 'all')
    rho_fine_IE
end

%% Evaluating points in s_patches
% refine one more time
rho_coarse_IE = rho_fine_IE;
gr_phi_coarse = gr_phi_fine;
curve_param_coarse = curve_param_fine;

rho_fine_IE = rho_coarse_IE + 1;
[gr_phi_fine, curve_param_fine] = refine_gr_phi(curve_param_rho1, rho_fine_IE, gr_phi_fft_rho1);

% excludes boundary
for eta_idx = M:-1:2
    for i = 1:curve_seq.n_curves
        s_patch = s_patches{i};
        [patch_X_s, patch_Y_s] = s_patch.xy_mesh;
        for xi_idx = 1:s_patch.n_xi
            x = patch_X_s(eta_idx, xi_idx); y = patch_Y_s(eta_idx, xi_idx);
            u_num_coarse_val = IE_curve_seq.u_num(x, y, curve_param_coarse, gr_phi_coarse);
            u_num_fine_val = IE_curve_seq.u_num(x, y, curve_param_fine, gr_phi_fine);
            
            while abs(u_num_coarse_val - u_num_fine_val) > int_eps
                rho_coarse_IE = rho_fine_IE
                curve_param_coarse = curve_param_fine;
                gr_phi_coarse = gr_phi_fine;

                rho_fine_IE = rho_coarse_IE + 1;
                [gr_phi_fine, curve_param_fine] = refine_gr_phi(curve_param_rho1, rho_fine_IE, gr_phi_fft_rho1);

                u_num_coarse_val = u_num_fine_val;
                u_num_fine_val = IE_curve_seq.u_num(x, y, curve_param_fine, gr_phi_fine);
            end
            s_patch.f_XY(eta_idx, xi_idx) = u_num_coarse_val;
        end
    end
end

% boundary value computation
for i = 1:curve_seq.n_curves
    s_patch = s_patches{i};
    [patch_X, patch_Y] = s_patch.xy_mesh;

    s_patch.f_XY(1, :) = u_G(patch_X(1, :), patch_Y(1, :));
end

for i = 1:curve_seq.n_curves
    s_patch = s_patches{i};
    in_patch = s_patch_msks{i};
    R_s_idxs = R.R_idxs(in_patch);

    [P_xi_s, P_eta_s] = R_xi_eta_inversion(R, s_patch, in_patch);

    for idx = 1:length(R_s_idxs)
        u_num_mat(R_s_idxs(idx)) = s_patch.locally_compute(P_xi_s(idx), P_eta_s(idx), M);
    end
end

err_msk = well_interior_msk;
for i = 1:curve_seq.n_curves
    err_msk = err_msk | s_patch_msks{i};
end



%% Evaluating points in corner patches using most refined gr_phi

curr = IE_curve_seq.first_curve;
for curve_idx = 1:IE_curve_seq.n_curves
    next = curr.next_curve;
    
    % C2-type corner case
    if curr.C2_corner
        % corner patches
        c_1_patch = c_1_patches{curve_idx};
        c_0_patch = c_0_patches{next.curve_idx};

        % computes windowed graded phis and the respective integration
        % intervals, these are the same for all of the points in the corner
        % patch corresponding to the same corner
        gr_phi_U_c_0 =  next.c_0_apply_POU(curve_param_fine, U, gr_phi_fine);
        U_c_0_idx_seg = [1, curve_param_fine.U_c_0_intervals(next.curve_idx, 2)-1];
        U_c_0_curve_seg = [next.curve_idx, next.curve_idx];

        gr_phi_U_c_0_C = apply_POU(gr_phi_fine, U, curve_param_fine.U_c_0_intervals(next.curve_idx, :));
        U_c_0_C_idx_seg = [curve_param_fine.U_c_0_intervals(next.curve_idx, 1), curve_param_fine.curve_n(curve_idx)];
        U_c_0_C_curve_seg = [next.curve_idx, curve_idx];

        gr_phi_U_c_1 = curr.c_1_apply_POU(curve_param_fine, U, gr_phi_fine);
        U_c_1_idx_seg = [curve_param_fine.U_c_1_intervals(curve_idx, 1), curve_param_fine.curve_n(curve_idx)];
        U_c_1_curve_seg = [curr.curve_idx, curr.curve_idx];

        gr_phi_U_c_1_C = apply_POU(gr_phi_fine, 1-U, curve_param_fine.U_c_1_intervals(curve_idx, :));
        U_c_1_C_idx_seg = [1, curve_param_fine.U_c_1_intervals(next.curve_idx, 2)-1];
        U_c_1_C_curve_seg = [next.curve_idx, curve_idx];

        gr_phi_U_c_0_1_C = apply_POU(gr_phi_U_c_1_C, U, curve_param_fine.U_c_0_intervals(next.curve_idx, :));
        U_c_0_1_C_idx_seg = [curve_param_fine.U_c_0_intervals(next.curve_idx, 1), curve_param_fine.U_c_1_intervals(curve_idx, 2)-1];
        U_c_0_1_C_curve_seg = [next.curve_idx, curve_idx];

        % Computing function values of cartesian points in c_0 and c_1 patches
        c_0_patch_msk = c_0_patch_msks{next.curve_idx};
        c_1_patch_msk = c_1_patch_msks{curve_idx};

        [~, h_eta_0] = c_0_patch.h_mesh;
        [~, h_eta_1] = c_1_patch.h_mesh;
        [P_xi_0, P_eta_0] = R_xi_eta_inversion(R, c_0_patch, c_0_patch_msk);
        [P_xi_1, P_eta_1] = R_xi_eta_inversion(R, c_1_patch, c_1_patch_msk);

        [R_c_0_intp_idxs, R_c_1_intp_idxs, R_c_0_1_intp_m_0_idxs, R_c_0_1_intp_m_1_idxs, R_no_intp_idxs] = construct_R_c_idxs(R, c_0_patch_msk, c_1_patch_msk, h_eta_0, h_eta_1, P_eta_0, P_eta_1);

        err_msk(R_no_intp_idxs) = true;
        for i = 1:length(R_no_intp_idxs)
            u_num_mat(R_no_intp_idxs(i)) = IE_curve_seq.u_num(R.R_X(R_no_intp_idxs(i)), R.R_Y(R_no_intp_idxs(i)), curve_param_fine, gr_phi_fine);
        end

        err_msk(R_c_0_intp_idxs(:, 1)) = true;
        for i = 1:size(R_c_0_intp_idxs, 1)
            u_num_mat(R_c_0_intp_idxs(i, 1)) = IE_curve_seq.interp_int_seg( ...
                P_xi_0(R_c_0_intp_idxs(i, 2)), ...
                P_eta_0(R_c_0_intp_idxs(i, 2)), ...
                c_0_patch, ...
                next, ...
                curve_param_fine, ...
                gr_phi_U_c_0, ...
                U_c_0_curve_seg, ...
                U_c_0_idx_seg, ...
                u_G, ...
                gr_phi_U_c_0_C, ...
                U_c_0_C_curve_seg, ...
                U_c_0_C_idx_seg ...
            );
            u_num_mat(R_c_0_intp_idxs(i, 1)) = u_num_mat(R_c_0_intp_idxs(i, 1)) + ...
                IE_curve_seq.int_segment_general( ...
                    curve_param_fine, ...
                    R.R_X(R_c_0_intp_idxs(i, 1)), ...
                    R.R_Y(R_c_0_intp_idxs(i, 1)), ...
                    gr_phi_U_c_0_C, ...
                    U_c_0_C_curve_seg, ...
                    U_c_0_C_idx_seg ...
                );
        end

        err_msk(R_c_1_intp_idxs(:, 1)) = true;
        for i = 1:size(R_c_1_intp_idxs, 1)
            u_num_mat(R_c_1_intp_idxs(i, 1)) = IE_curve_seq.interp_int_seg( ...
                P_xi_1(R_c_1_intp_idxs(i, 2)), ...
                P_eta_1(R_c_1_intp_idxs(i, 2)), ...
                c_1_patch, ...
                curr, ...
                curve_param_fine, ...
                gr_phi_U_c_1, ...
                U_c_1_curve_seg, ...
                U_c_1_idx_seg, ...
                u_G, ...
                gr_phi_U_c_1_C, ...
                U_c_1_C_curve_seg, ...
                U_c_1_C_idx_seg ...
            );
            u_num_mat(R_c_1_intp_idxs(i, 1)) = u_num_mat(R_c_1_intp_idxs(i, 1)) + ...
                IE_curve_seq.int_segment_general( ...
                    curve_param_fine, ...
                    R.R_X(R_c_1_intp_idxs(i, 1)), ...
                    R.R_Y(R_c_1_intp_idxs(i, 1)), ...
                    gr_phi_U_c_1_C, ...
                    U_c_1_C_curve_seg, ...
                    U_c_1_C_idx_seg ...
                );
        end

        for i = 1:size(R_c_0_1_intp_m_0_idxs, 1)
            u_num_mat(R_c_0_1_intp_m_0_idxs(i, 1)) = IE_curve_seq.interp_int_seg( ...
                P_xi_0(R_c_0_1_intp_m_0_idxs(i, 2)), ...
                P_eta_0(R_c_0_1_intp_m_0_idxs(i, 2)), ...
                c_0_patch, ...
                next, ...
                curve_param_fine, ...
                gr_phi_U_c_0, ...
                U_c_0_curve_seg, ...
                U_c_0_idx_seg, ...
                u_G, ...
                gr_phi_U_c_0_C, ...
                U_c_0_C_curve_seg, ...
                U_c_0_C_idx_seg ...
            );
            u_num_mat(R_c_0_1_intp_m_1_idxs(i, 1)) = u_num_mat(R_c_0_1_intp_m_1_idxs(i, 1)) +...
                IE_curve_seq.interp_int_seg( ...
                P_xi_1(R_c_0_1_intp_m_1_idxs(i, 2)), ...
                P_eta_1(R_c_0_1_intp_m_1_idxs(i, 2)), ...
                c_1_patch, ...
                curr, ...
                curve_param_fine, ...
                gr_phi_U_c_1, ...
                U_c_1_curve_seg, ...
                U_c_1_idx_seg, ...
                u_G, ...
                gr_phi_U_c_1_C, ...
                U_c_1_C_curve_seg, ...
                U_c_1_C_idx_seg ...
            );
            u_num_mat(R_c_0_1_intp_m_0_idxs(i, 1)) = u_num_mat(R_c_0_1_intp_m_0_idxs(i, 1)) + ...
                IE_curve_seq.int_segment_general( ...
                    curve_param_fine, ...
                    R.R_X(R_c_0_1_intp_m_0_idxs(i, 1)), ...
                    R.R_Y(R_c_0_1_intp_m_0_idxs(i, 1)), ...
                    gr_phi_U_c_0_1_C, ...
                    U_c_0_1_C_curve_seg, ...
                    U_c_0_1_C_idx_seg ...
                );
        end
    else
    end
    curr = curr.next_curve;
end
u_num_mat(~R.in_interior) = nan;
max(u_num_mat(:) - u_G(R.R_X(:), R.R_Y(:)))
end

function [gr_phi_rho, curve_param_rho] = refine_gr_phi(curve_param_rho1, rho_IE, gr_phi_fft_rho1)
    curve_param_rho = curve_param_obj(curve_param_rho1.curve_n * rho_IE);
    
    padded_fft_coeffs = [zeros(ceil((curve_param_rho1.n_total*rho_IE-curve_param_rho1.n_total)/2), 1); gr_phi_fft_rho1; zeros(floor((curve_param_rho1.n_total*rho_IE-curve_param_rho1.n_total)/2), 1)];
    gr_phi_rho = rho_IE*curve_param_rho1.n_total*real(ifft(ifftshift(padded_fft_coeffs)));
end

function [well_interior_msk, s_patch_msks, c_0_patch_msks, c_1_patch_msks] = gen_R_msks(R, n_r, s_patches, c_0_patches, c_1_patches)
    well_interior_msk = R.in_interior;
    s_patch_msks = cell(size(s_patches));
    c_0_patch_msks = cell(size(c_0_patches));
    c_1_patch_msks = cell(size(c_1_patches));
    
    for i = 1:length(s_patch_msks)
        %S_patch
        s_patch = s_patches{i};
        [bound_X, bound_Y] = s_patch.boundary_mesh_xy(n_r, false);
        in_patch = inpolygon_mesh(R.R_X, R.R_Y, bound_X, bound_Y) & R.in_interior;
        s_patch_msks{i} = in_patch;
        well_interior_msk = well_interior_msk & ~in_patch;

        c_0_patch = c_0_patches{i};
        [bound_X, bound_Y] = c_0_patch.boundary_mesh_xy(n_r, false);
        in_patch = inpolygon_mesh(R.R_X, R.R_Y, bound_X, bound_Y) & R.in_interior;
        c_0_patch_msks{i} = in_patch;
        well_interior_msk = well_interior_msk & ~in_patch;

        c_1_patch = c_1_patches{i};
        [bound_X, bound_Y] = c_1_patch.boundary_mesh_xy(n_r, false);
        in_patch = inpolygon_mesh(R.R_X, R.R_Y, bound_X, bound_Y) & R.in_interior;
        c_1_patch_msks{i} = in_patch;
        well_interior_msk = well_interior_msk & ~in_patch;
    end
end

function U = construct_POU(n_POU)
    w_1D = @(x) erfc(6*(-2*x+1))/2;
    w_1 = @(s) w_1D(s);
    w_2 = @(s) w_1D(-s+1);
    U_mesh = linspace(0, 1, n_POU)';
    
    U = w_1(U_mesh)./(w_1(U_mesh)+w_2(U_mesh));
end

function gr_phi_POU = apply_POU(gr_phi, U, POU_interval)
gr_phi_POU = gr_phi; gr_phi_POU(POU_interval(1):POU_interval(2)) = gr_phi_POU(POU_interval(1):POU_interval(2)) .* U;
end

function [R_c_0_intp_idxs, R_c_1_intp_idxs, R_c_0_1_intp_m_0_idxs, R_c_0_1_intp_m_1_idxs, R_no_intp_idxs] = construct_R_c_idxs(R, c_0_patch_msk, c_1_patch_msk, h_eta_0, h_eta_1, P_eta_0, P_eta_1)
    c_0_all_idxs = (1:sum(c_0_patch_msk, 'all'))';
    c_1_all_idxs = (1:sum(c_1_patch_msk, 'all'))';
    R_c_0_all_idxs = R.R_idxs(c_0_patch_msk);
    R_c_1_all_idxs = R.R_idxs(c_1_patch_msk);
    
    intp_c_0_all_msk = false(size(R.f_R)); intp_c_0_all_msk(R_c_0_all_idxs(P_eta_0 < h_eta_0)) = true;
    intp_c_1_all_msk = false(size(R.f_R)); intp_c_1_all_msk(R_c_1_all_idxs(P_eta_1 < h_eta_1)) = true;
    R_no_intp_idxs = R.R_idxs((c_0_patch_msk | c_1_patch_msk) & ~(intp_c_0_all_msk | intp_c_1_all_msk));
    
    intp_c_0_msk = intp_c_0_all_msk & ~intp_c_1_all_msk;
    intp_c_1_msk = intp_c_1_all_msk & ~intp_c_0_all_msk;
    intp_c_0_1_msk = intp_c_0_all_msk & intp_c_1_all_msk;
    R_c_0_intp_idxs = [R_c_0_all_idxs(intp_c_0_msk(R_c_0_all_idxs)), c_0_all_idxs(intp_c_0_msk(R_c_0_all_idxs))];
    R_c_1_intp_idxs = [R_c_1_all_idxs(intp_c_1_msk(R_c_1_all_idxs)), c_1_all_idxs(intp_c_1_msk(R_c_1_all_idxs))];
    R_c_0_1_intp_m_0_idxs = [R_c_0_all_idxs(intp_c_0_1_msk(R_c_0_all_idxs)), c_0_all_idxs(intp_c_0_1_msk(R_c_0_all_idxs))];
    R_c_0_1_intp_m_1_idxs = [R_c_1_all_idxs(intp_c_0_1_msk(R_c_1_all_idxs)), c_1_all_idxs(intp_c_0_1_msk(R_c_1_all_idxs))];
end

