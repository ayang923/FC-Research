clc; clear; close all;
% 2DFC on boomerang domain
%% Setting Parameters
f = @(x, y) 4 + (1 + x.^2 + y.^2).*(sin(2.5*pi*x - 0.5) + cos(2*pi*y - 0.5));

n_C1_bound = 80;

C = 27;
d = 5;
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
load(['boomerang_data/patches_nC1bound', num2str(n_C1_bound), '_d', num2str(d)])
load(['boomerang_data/interior_mesh_nC1bound', num2str(n_C1_bound)])

%% Computing window functions for C1 Patch
[C1_B_norm, C1_B_c_norm, xi_norm, eta_norm] = C1_patch.compute_phi_normalization(window_patch_xi, window_patch_eta);

%% FC in Parameter Space
[C1_fcont_patch_xi, C1_fcont_patch_eta_refined, C1_fcont_patch_eta_unrefined] = C1_patch.FC(C, n_r, d, A, Q, M, C1_B_norm, C1_B_c_norm);
window_fcont_patch_xi = window_patch_xi.FC(C, n_r, d, A, Q, xi_norm);
window_fcont_patch_eta = window_patch_eta.FC(C, n_r, d, A, Q, eta_norm);
S_fcont_patch = S_patch.FC(C, n_r, d, A, Q, nan);

fcont_patches = {C1_fcont_patch_xi, C1_fcont_patch_eta_refined, C1_fcont_patch_eta_unrefined, window_fcont_patch_xi, window_fcont_patch_eta, S_fcont_patch};
%% Interpolation onto Cartesian Mesh
% computes bounds of R
R_x_bounds = [S_fcont_patch.x_min-h_R, S_fcont_patch.x_max+h_R];
R_y_bounds = [S_fcont_patch.y_min-h_R, S_fcont_patch.y_max+h_R];

% Constructs continued cartesian mesh object
R_FC = R.extend_R(R_x_bounds(1), R_x_bounds(2), R_y_bounds(1), R_y_bounds(2), true);

clear 'R'

%Fills interior with exact values and interpolates patch values onto grid
for patch = fcont_patches
    R_FC.interpolate_patch(patch{1}, d+3);
end

%% Plots
% figure;
% s = surf(R_FC.R_X, R_FC.R_Y, R_FC.f_R);
% s.EdgeColor = 'none';
% 
% figure;
% scatter(R_FC.R_X(:), R_FC.R_Y(:), 100, R_FC.f_R(:), 'filled');
% colorbar;
% colormap(jet);
% hold on;
% plot(R_FC.boundary_X, R_FC.boundary_Y)
% 
% figure;
% scatter3(R_FC.R_X(:), R_FC.R_Y(:), R_FC.f_R(:));

%% FFT and Error Calculation
R_FC.compute_fc_coeffs()
[X_err, Y_err, f_err, interior_idxs] = R_FC.ifft_interpolation(R_FC.h * 0.5, true);

% load finer interior data
% load(['boomerang_data/interior_mesh_nC1bound', num2str(n_C1_bound*2)])
% R_exact = R.extend_R(X_err(1, 1), X_err(1, end), Y_err(1, 1), Y_err(end, 1), false);
% clear 'R';

err = abs(f(X_err(interior_idxs), Y_err(interior_idxs)) - f_err(interior_idxs));
max(err, [], 'all') %l_infinity error

figure;
scatter(X_err(interior_idxs), Y_err(interior_idxs), 100, err, 'filled');

% Add a color bar
colorbar;

% Set colormap
colormap(jet);
