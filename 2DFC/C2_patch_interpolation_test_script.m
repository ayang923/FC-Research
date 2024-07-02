clc; clear; close all;

n_lst = 20 * 2.^(4:5);
h_lst = 0.029 * (1/2).^(4:5);

load("C2_patches.mat")
% 2D
C = 27;
d = 4;

load(['FC_data/A_d',num2str(d),'_C', num2str(27), '.mat']);
load(['FC_data/Q_d',num2str(d),'_C', num2str(27), '.mat']);
A = double(A);
Q = double(Q);

% 1D-FC on each cross section using finest mesh for error

% Computing windowed function on finest mesh
finest_C2_patch = C2_patches{length(n_lst)};
[XI_err, ETA_err] = finest_C2_patch.xi_eta_mesh();

finest_C2_patch.f_XY = finest_C2_patch.phi(XI_err, ETA_err) ./ C2_norms{length(n_lst)};
finest_C2_fcont_patch = finest_C2_patch.FC(C, d, A, Q, nan);
test_patch_i = 1;
% FC and FFT
C2_patch = C2_patches{test_patch_i};
C2_norm = C2_norms{test_patch_i};

[XI, ETA] = C2_patch.xi_eta_mesh();
% C2_patch.f_XY = C2_patch.f_XY .* C2_patch.phi(XI, ETA);

C2_patch.f_XY = C2_patch.phi(XI, ETA)./C2_norm;
C2_fcont_patch = C2_patch.FC(C, d, A, Q, nan);

[XI_fcont, ETA_fcont] = C2_fcont_patch.xi_eta_mesh();

[bound_C2_X, bound_C2_Y] = C2_patch.boundary_mesh_xy();

R_x_bounds = [C2_fcont_patch.x_min, C2_fcont_patch.x_max];
R_y_bounds = [C2_fcont_patch.y_min, C2_fcont_patch.y_max];

R = R_cartesian_mesh_obj(R_x_bounds(1), R_x_bounds(2), R_y_bounds(1), R_y_bounds(2), h_lst(1)./2, bound_C2_X, bound_C2_Y);
R_err = R_cartesian_mesh_obj(R_x_bounds(1), R_x_bounds(2), R_y_bounds(1), R_y_bounds(2), h_lst(1)./2, bound_C2_X, bound_C2_Y);
figure;
scatter(R.R_X(:), R.R_Y(:))
hold on;
scatter(R.R_X(R.in_interior), R.R_Y(R.in_interior))

R.interpolate_patch(C2_fcont_patch, d+3, false);
R_err.interpolate_patch(finest_C2_patch, d+3, false);