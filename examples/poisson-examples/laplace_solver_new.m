function [u_num_mat] = laplace_solver_new(R, curve_seq, u_G, cfac, p, M, int_eps, eps_xi_eta, eps_xy, n_r)
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

[s_patches, c_0_patches, c_1_patches] = IE_curve_seq.construct_interior_patches(curve_param_rho1, R.h, M, eps_xi_eta, eps_xy);
[well_interior_msk, s_patch_msks, c_0_patch_msks, c_1_patch_msks] = gen_R_msks(R, n_r, s_patches, c_0_patches, c_1_patches);

test_curve = IE_curve_seq.first_curve;
[gr_phi_U] = test_curve.c_1_apply_POU(curve_param_rho1, U, gr_phi_rho1);
u_num_mat = nan;

x = 1e-3;
y = 0.6e-3;

rho_fine_IE = 4;
[gr_phi_fine, curve_param_fine] = refine_gr_phi(curve_param_rho1, rho_fine_IE, gr_phi_fft_rho1);

i = 1;
curr = IE_curve_seq.first_curve;

[xi, ~, ~] = c_1_patches{i}.inverse_M_p(x, y, nan);
eta = 1e-6;

x = M_p_x(c_1_patches{i}, xi, eta);
y = M_p_y(c_1_patches{i}, xi, eta);

eta_mesh = c_1_patches{i}.eta_mesh*1/32;
interpol_nodes = c_1_patches{i}.M_p(xi, eta_mesh);

interpol_vals = zeros(size(eta_mesh));

test_curve = IE_curve_seq.first_curve;
gr_phi_fine_U = test_curve.c_1_apply_POU(curve_param_fine, U, gr_phi_fine);
seg_interval = [curve_param_fine.U_c_1_intervals(test_curve.curve_idx, 1), curve_param_fine.curve_n(test_curve.curve_idx)];
for i = 1:length(interpol_vals)
    interpol_vals(i) = IE_curve_seq.int_segment_general(curve_param_fine, interpol_nodes(i, 1), interpol_nodes(i, 2), gr_phi_fine_U, [1, 1], seg_interval);
end

gr_phi_fine_U_complement = apply_POU(gr_phi_fine, 1-U, curve_param_fine.U_c_1_intervals(test_curve.curve_idx, :));

seg_interval_complement = [1, curve_param_fine.U_c_1_intervals(test_curve.curve_idx, 2)-1];
interpol_vals(1) = u_G(interpol_nodes(1, 1), interpol_nodes(1, 2)) - IE_curve_seq.int_segment_boundary(curve_param_fine, xi, test_curve, gr_phi_fine_U_complement, [2, 1], seg_interval_complement);

int_contr_hat = barylag([eta_mesh, interpol_vals], eta);
% int_contr_hat = IE_curve_seq.int_segment_general(curve_param_fine, x, y, gr_phi_fine_U, [1, 1], seg_interval)
int_contr = IE_curve_seq.int_segment_general(curve_param_fine, x, y, gr_phi_fine_U, [1, 1], seg_interval);

figure;
fplot(@(X) arrayfun(@(x) barylag([eta_mesh, interpol_vals], x), X), [eta_mesh(1), eta_mesh(end)])
hold on;
plot(eta_mesh, interpol_vals, 'x')
plot(eta, int_contr, 'x');

int_seg_fun = @(eta_tilde) IE_curve_seq.int_segment_general(curve_param_fine, M_p_x(c_1_patches{1}, xi, eta_tilde), M_p_y(c_1_patches{1}, xi, eta_tilde), gr_phi_fine_U, [1, 1], seg_interval);

% fplot(@(ETA) arrayfun(@(eta) int_seg_fun(eta), ETA), [eta_mesh(1), eta_mesh(end)])

est_hat = int_contr_hat + IE_curve_seq.int_segment_general(curve_param_fine, x, y, gr_phi_fine_U_complement, [1, 1], seg_interval_complement);
est = int_contr + IE_curve_seq.int_segment_general(curve_param_fine, x, y, gr_phi_fine_U_complement, [1, 1],  seg_interval_complement);

est_hat - (x.^2 - y.^2)
est - (x.^2-y.^2)

% %% Fills in points that are well within the domain
% u_num_mat = zeros(size(R.f_R));
% u_num_mat_fine = zeros(size(R.f_R));
% 
% for idx = R.R_idxs(well_interior_msk)'
%     u_num_mat(idx) = IE_curve_seq.u_num(R.R_X(idx), R.R_Y(idx), curve_param_rho1, gr_phi_rho1);
% end
% 
% to_update = well_interior_msk;
% rho_fine_IE = 2;
% 
% while sum(to_update, 'all') > 0
%     [gr_phi_fine, curve_param_fine] = refine_gr_phi(curve_param_rho1, rho_fine_IE, gr_phi_fft_rho1);
% 
%     for idx = R.R_idxs(to_update)'
%         u_num_mat_fine(idx) = IE_curve_seq.u_num(R.R_X(idx), R.R_Y(idx), curve_param_fine, gr_phi_fine);
%     end
% 
%     to_update = to_update & abs(u_num_mat_fine - u_num_mat) > int_eps;
%     u_num_mat(to_update) = u_num_mat_fine(to_update);
%     rho_fine_IE = rho_fine_IE + 1;
% 
%     sum(to_update, 'all')
%     rho_fine_IE
% end

%% Evaluating points in s_patches
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

function M_p_x = M_p_x(patch, xi, eta)
    M_p_xy = patch.M_p(xi, eta);
    M_p_x = M_p_xy(1);
end
function M_p_y = M_p_y(patch, xi, eta)
    M_p_xy = patch.M_p(xi, eta);
    M_p_y = M_p_xy(2);
end
