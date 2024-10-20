clc; clear; close all;
% 2DFC on teardrop domain
%% Setting Parameters
f = @(x, y) 4 + (1 + x.^2 + y.^2).*(sin(2.5*pi*x - 0.5) + cos(2*pi*y - 0.5));
l_theta = @(theta) [2*sin(theta/2), -sin(theta)];

scale_factor = 2;

n_C2 = 40*scale_factor;

C = 27;
d = 7;
n_r = 6;

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
load(['teardrop_data/patches_nC2', num2str(n_C2), '_d', num2str(d)])
load(['teardrop_data/interior_mesh_nC2', num2str(n_C2)])

%% Computing window functions for C2 Patch
[C2_norm, xi_norm, eta_norm] = C2_patch.compute_phi_normalization(window_patch_xi, window_patch_eta);

%% FC in Parameter Space
[C2_fcont_patch_xi, C2_fcont_patch_eta, C2_fcont_patch_corner] = C2_patch.FC(C, n_r, d, A, Q, C2_norm);
window_fcont_patch_xi = window_patch_xi.FC(C, n_r, d, A, Q, xi_norm);
window_fcont_patch_eta = window_patch_eta.FC(C, n_r, d, A, Q, eta_norm);
S_fcont_patch = S_patch.FC(C, n_r, d, A, Q, nan);

fcont_patches = {C2_fcont_patch_xi, C2_fcont_patch_eta, C2_fcont_patch_corner, window_fcont_patch_xi, window_fcont_patch_eta, S_fcont_patch};

%% Interpolation onto Cartesian Mesh
% computes bounds of R
R_x_bounds = [C2_fcont_patch_corner.x_min-R.h*2, S_fcont_patch.x_max+R.h*2];
R_y_bounds = [S_fcont_patch.y_min-R.h*2, S_fcont_patch.y_max+R.h*2];

% Constructs cartesian mesh object
R_FC = R.extend_R(R_x_bounds(1), R_x_bounds(2), R_y_bounds(1), R_y_bounds(2), true);

clear 'R';

%Fills interior with exact values and interpolates patch values onto grid
for patch = fcont_patches
    R_FC.interpolate_patch(patch{1}, d+3);
end
R_FC.fill_interior(f);

%% Plots
figure;
s = surf(R_FC.R_X, R_FC.R_Y, R_FC.f_R);
s.EdgeColor = 'none';

figure;
scatter3(R_FC.R_X(:), R_FC.R_Y(:), R_FC.f_R(:));

%% FFT and Error Calculation
R_FC.compute_fc_coeffs()
[X_err, Y_err, f_err, err_interior_idx] = R_FC.ifft_interpolation(R_FC.h * 0.5, true);

max(abs((f(X_err(err_interior_idx), Y_err(err_interior_idx)) - f_err(err_interior_idx)))) %l_infinity error

err = abs((f(X_err(err_interior_idx), Y_err(err_interior_idx)) - f_err(err_interior_idx)));

figure;
scatter(X_err(err_interior_idx), Y_err(err_interior_idx), 100, err, 'filled');

% Add a color bar
colorbar;

% Set colormap
colormap(jet);
