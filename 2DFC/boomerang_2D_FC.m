clc; clear; close all;
% 2DFC on boomerang domain
%% Setting Parameters
f = @(x, y) 4 + (1 + x.^2 + y.^2).*(sin(2.5*pi*x - 0.5) + cos(2*pi*y - 0.5));

n_C1_bound = 40;

C = 27;
d = 4;
n_r = 6;

M = d+3;

%% Loading FC Data
if(exist(['FC_data/A_d',num2str(d),'_C', num2str(C), '_r', num2str(n_r), '.mat']) == 0 | ...
   exist(['FC_data/Q_d',num2str(d),'_C', num2str(C),  '_r', num2str(n_r), '.mat']) == 0)
    disp('FC data not found. Generating FC operators... \n');
    generate_bdry_continuations(d, C, C, 12, 20, 4, ...
        256, n_r);
end

load(['FC_data/A_d',num2str(d),'_C', num2str(C), '_r', num2str(n_r), '.mat']);
load(['FC_data/Q_d',num2str(d),'_C', num2str(C),  '_r', num2str(n_r), '.mat']);

A = double(A);
Q = double(Q);

%% Loading Precomputed Patches and Interior Mesh
load(['boomerang_data/boomerang_data_nC1bound', num2str(n_C1_bound)])

%% Computing window functions for C1 Patch
[C1_B_norm, C1_B_c_norm, xi_norm, eta_norm] = C1_patch.compute_phi_normalization(window_patch_xi, window_patch_eta);

figure;
[XI, ETA] = C1_patch.B.xi_eta_mesh;
surf(XI, ETA, C1_patch.B.f_XY.*C1_patch.B.phi(XI, ETA)./C1_B_norm)

%% FC in Parameter Space
[C1_fcont_patch_xi, C1_fcont_patch_eta_refined, C1_fcont_patch_eta_unrefined] = C1_patch.FC(C, n_r, d, A, Q, M, C1_B_norm, C1_B_c_norm);
window_fcont_patch_xi = window_patch_xi.FC(C, n_r, d, A, Q, xi_norm);
window_fcont_patch_eta = window_patch_eta.FC(C, n_r, d, A, Q, eta_norm);
S_fcont_patch = S_patch.FC(C, n_r, d, A, Q, nan);

fcont_patches = {C1_fcont_patch_xi, C1_fcont_patch_eta_refined, C1_fcont_patch_eta_unrefined, window_fcont_patch_xi, window_fcont_patch_eta, S_fcont_patch};
% fcont_patches = {C1_fcont_patch_xi, C1_fcont_patch_eta_refined, C1_fcont_patch_eta_unrefined};

%% Interpolation onto Cartesian Mesh
% computes bounds of R
R_x_bounds = [S_fcont_patch.x_min-h_R, S_fcont_patch.x_max+h_R];
R_y_bounds = [S_fcont_patch.y_min-h_R, S_fcont_patch.y_max+h_R];

% Constructs continued cartesian mesh object
R = R.extend_R(R_x_bounds(1), R_x_bounds(2), R_y_bounds(1), R_y_bounds(2));

%Fills interior with exact values and interpolates patch values onto grid
for patch = fcont_patches
    R.interpolate_patch(patch{1}, d+3, false);
end
%% Plots
figure;
s = surf(R.R_X, R.R_Y, R.f_R);
s.EdgeColor = 'none';

figure;
scatter(R.R_X(:), R.R_Y(:), 100, R.f_R(:), 'filled');
colorbar;
colormap(jet);
hold on;
plot(R.boundary_X, R.boundary_Y)

figure;
scatter3(R.R_X(:), R.R_Y(:), R.f_R(:));

%% FFT and Error Calculation
R.compute_fc_coeffs()
[X_err, Y_err, f_err, err_interior_idx] = R.ifft_interpolation(R.h * 0.5);

max(abs((f(X_err(err_interior_idx), Y_err(err_interior_idx)) - f_err(err_interior_idx)))) %l_infinity error

err = abs((f(X_err(err_interior_idx), Y_err(err_interior_idx)) - f_err(err_interior_idx)));

figure;
scatter(X_err(err_interior_idx), Y_err(err_interior_idx), 100, err, 'filled');

% Add a color bar
colorbar;

% Set colormap
colormap(jet);
