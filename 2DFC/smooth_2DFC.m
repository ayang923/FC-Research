clc; clear; close all;

%% Setting Parameters
f = @(x, y) -(x.^6+y.^6).*sin(10*pi*x).*sin(10*pi*y);
l_theta = @(theta) [cos(theta)+0.35.*cos(2.*theta)-0.35, 0.7.*sin(theta)];


scale_factor = 4;

h = 0.01 / scale_factor;
n_S = round(pi/h);

C = 27;
d = 4;
n_r = 6;

M = d+3;

%% Creating Patches
S_patch_1 = construct_S_patch(f, 0, pi, h, n_S, d);
S_patch_2 = construct_S_patch(f, pi, 2*pi, h, n_S, d);

figure;
[X, Y] = S_patch_1.xy_mesh();
scatter(X(:), Y(:))
hold on;
[X, Y] = S_patch_2.xy_mesh();
scatter(X(:), Y(:))
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

%% FC
S_fcont_patch_1 = S_patch_1.FC(C, n_r, d, A, Q, nan);
S_fcont_patch_2 = S_patch_2.FC(C, n_r, d, A, Q, nan);

fcont_patches = {S_fcont_patch_1, S_fcont_patch_2};

%% Interpolation onto Cartesian Mesh
% computes bounds of R
R_x_bounds = [S_fcont_patch_1.x_min-0.01, S_fcont_patch_1.x_max+0.01];
R_y_bounds = [S_fcont_patch_2.y_min-0.01, S_fcont_patch_1.y_max+0.01];

% Computes boundary_XY values of domain 
boundary_XY = l_theta(transpose(linspace(0, 2*pi, 200*scale_factor)));

% Constructs cartesian mesh object
R = R_cartesian_mesh_obj(R_x_bounds(1), R_x_bounds(2), R_y_bounds(1), R_y_bounds(2), h, boundary_XY(:, 1), boundary_XY(:, 2));

%Fills interior with exact values and interpolates patch values onto grid
for patch = fcont_patches
    R.interpolate_patch(patch{1}, M, false, f)
end
R.fill_interior(f);

figure;
s = surf(R.R_X, R.R_Y, R.f_R);
s.EdgeColor = 'none';

figure;
scatter3(R.R_X(R.interior_idxs), R.R_Y(R.interior_idxs), R.f_R(R.interior_idxs));

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

function S_patch = construct_S_patch(f, theta_A, theta_B, h, n_xi, n_eta)
    l_theta = @(theta) [cos(theta)+0.35.*cos(2.*theta)-0.35, 0.7.*sin(theta)];
    
    l_A = @(xi) (theta_B - theta_A).*xi + theta_A;
    nu =  @(theta) [-0.7.*cos(theta), -sin(theta)-0.7.*sin(2.*theta)] ./ sqrt((0.7.*cos(theta)).^2 + (sin(theta)+0.7.*sin(2.*theta)).^2);
    
    M_p_general = @(xi, eta, H) l_theta(l_A(xi)) + eta.*H.*nu(l_A(xi));
    % H is a function of h and n_eta
    
    theta_diff = theta_B - theta_A;
    
    dM_pdxi = @(xi, eta, H) [
      -theta_diff .* sin(theta_A + theta_diff .* xi) + ...
      (0.7 .* eta .* H .* theta_diff .* sin(theta_A + theta_diff .* xi)) ./ ...
      sqrt(0.49 .* cos(theta_A + theta_diff .* xi)^2 + ...
      (sin(theta_A + theta_diff .* xi) + 0.7 .* sin(2 .* (theta_A + theta_diff .* xi)))^2) - ...
      0.7 .* theta_diff .* sin(2 .* (theta_A + theta_diff .* xi)) + ...
      (0.35 .* eta .* H .* cos(theta_A + theta_diff .* xi) .* ...
      (2 .* theta_diff .* (cos(theta_A + theta_diff .* xi) + 1.4 .* cos(2 .* (theta_A + theta_diff .* xi))) .* ...
      (sin(theta_A + theta_diff .* xi) + 0.7 .* sin(2 .* (theta_A + theta_diff .* xi))) - ...
      0.49 .* theta_diff .* sin(2 .* (theta_A + theta_diff .* xi)))) ./ ...
      (0.49 .* cos(theta_A + theta_diff .* xi)^2 + ...
      (sin(theta_A + theta_diff .* xi) + 0.7 .* sin(2 .* (theta_A + theta_diff .* xi)))^2)^(3./2); ...

      0.7 .* theta_diff .* cos(theta_A + theta_diff .* xi) - ...
      (eta .* H .* theta_diff .* (cos(theta_A + theta_diff .* xi) + ...
      1.4 .* cos(2 .* (theta_A + theta_diff .* xi)))) ./ ...
      sqrt(0.49 .* cos(theta_A + theta_diff .* xi)^2 + ...
      (sin(theta_A + theta_diff .* xi) + 0.7 .* sin(2 .* (theta_A + theta_diff .* xi)))^2) - ...
      (eta .* H .* (-sin(theta_A + theta_diff .* xi) - 0.7 .* sin(2 .* (theta_A + theta_diff .* xi))) .* ...
      (2 .* theta_diff .* (cos(theta_A + theta_diff .* xi) + 1.4 .* cos(2 .* (theta_A + theta_diff .* xi))) .* ...
      (sin(theta_A + theta_diff .* xi) + 0.7 .* sin(2 .* (theta_A + theta_diff .* xi))) - ...
      0.49 .* theta_diff .* sin(2 .* (theta_A + theta_diff .* xi)))) ./ ...
      (2 .* (0.49 .* cos(theta_A + theta_diff .* xi)^2 + ...
      (sin(theta_A + theta_diff .* xi) + 0.7 .* sin(2 .* (theta_A + theta_diff .* xi)))^2)^(3./2))
    ];

    dM_pdeta = @(xi, eta, H) [
      -(0.7 .* H .* cos(theta_A + theta_diff .* xi) ./ ...
        sqrt(0.49 .* cos(theta_A + theta_diff .* xi)^2 + ...
             (sin(theta_A + theta_diff .* xi) + 0.7 .* sin(2 .* (theta_A + theta_diff .* xi)))^2)); ...
      -(H .* (sin(theta_A + theta_diff .* xi) + 0.7 .* sin(2 .* (theta_A + theta_diff .* xi))) ./ ...
        sqrt(0.49 .* cos(theta_A + theta_diff .* xi)^2 + ...
             (sin(theta_A + theta_diff .* xi) + 0.7 .* sin(2 .* (theta_A + theta_diff .* xi)))^2))
    ];

    J = @(v, H) [dM_pdxi(v(1), v(2), H), dM_pdeta(v(1), v(2), H)];
    
    xi_mesh = linspace(0, 1, n_xi+1);
    eta_mesh = linspace(0, 1, n_eta+1);
    
    [XI, ETA] = meshgrid(xi_mesh, eta_mesh);
    
    XY = M_p_general(XI(:), ETA(:), h.*n_eta);
    
    f_XY = reshape(f(XY(:, 1), XY(:, 2)), [n_eta+1, n_xi+1]);
    
    S_patch = S_patch_obj(M_p_general, J, n_xi, n_eta, 0, 1, 0, 1, f_XY, h, nan);
end