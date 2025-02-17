clc; clear; close all;
% 2DFC on teardrop domain
%% Setting Parameters
f = @(x, y) 4 + (1 + x.^2 + y.^2).*(sin(2.5*pi*x - 0.5) + cos(2*pi*y - 0.5));
l_theta = @(theta) [2*sin(theta/2), -sin(theta)];

scale_factor = 4;

h_R = 0.01 / scale_factor;

n_C2 = round(0.8/(h_R*2));
n_S_h = round(0.4/h_R);
h_S = h_R;

n_S_w = round((2*pi-1)/h_R);

C = 27;
d = 4;
n_r = 8;

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

%% Constructing Boundary Patches

C2_patch = construct_C2_patch(f, 0.4, 2*pi-0.4, 0, n_C2, n_C2); % data associated with patch

window_patch_xi = construct_S_patch(f, 0.1, 0.5, h_S, n_S_h, d+10);
window_patch_eta = construct_S_patch(f, 2*pi-0.1, 2*pi-0.5, h_S, n_S_h, d+10);
S_patch = construct_S_patch(f, 0.5, 2*pi-0.5, h_S, n_S_w, d+10);

figure;
[X, Y] = S_patch.xy_mesh();
scatter(X(:), Y(:))
hold on;
[X, Y] = C2_patch.xy_mesh();
scatter(X(:), Y(:))
[X, Y] = window_patch_xi.xy_mesh();
scatter(X(:), Y(:))
[X, Y] = window_patch_eta.xy_mesh();
scatter(X(:), Y(:))

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
R_x_bounds = [C2_fcont_patch_corner.x_min, S_fcont_patch.x_max];
R_y_bounds = [S_fcont_patch.y_min, S_fcont_patch.y_max];

% Computes boundary_XY values of domain 
boundary_XY = l_theta(transpose(linspace(0, 2*pi, 200*scale_factor)));

% Constructs cartesian mesh object
R = R_cartesian_mesh_obj(R_x_bounds(1), R_x_bounds(2), R_y_bounds(1), R_y_bounds(2), h_R, boundary_XY(:, 1), boundary_XY(:, 2));

%Fills interior with exact values and interpolates patch values onto grid
for patch = fcont_patches
    R.interpolate_patch(patch{1}, d+3, false, f);
end
R.fill_interior(f);

figure;
s = surf(R.R_X, R.R_Y, R.f_R);
s.EdgeColor = 'none';

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
%% Problem specific patch construction
function C2_patch = construct_C2_patch(f, theta_A, theta_B, theta_C, n_xi, n_eta)
    l_theta = @(theta) [2*sin(theta/2), -sin(theta)];
    
    l_A = @(xi) xi*(theta_A - theta_C) + theta_C;
    l_B = @(eta) eta*(theta_B - theta_C - 2*pi) + theta_C + 2*pi; % account for wrap around
    
    M_p = @(xi, eta) l_theta(l_A(xi)) + l_theta(l_B(eta)) - l_theta(theta_C);
    J = @(v) [theta_A*cos(v(1)*theta_A/2) (theta_B-2*pi)*cos((v(2)*(theta_B-2*pi)+2*pi)/2); -theta_A*cos(v(1)*theta_A) -(theta_B-2*pi)*cos(v(2)*(theta_B-2*pi)+2*pi)];
    
    xi_mesh = linspace(0, 1, n_xi+1);
    eta_mesh = linspace(0, 1, n_eta+1);
    
    [XI, ETA] = meshgrid(xi_mesh, eta_mesh);
    XY = M_p(XI(:), ETA(:));
    
    f_XY = reshape(f(XY(:, 1), XY(:, 2)), [n_eta+1, n_xi+1]);
    
    C2_patch = C2_patch_obj(M_p, J, n_xi, n_eta, 0, 1, 0, 1, f_XY, nan); % data associated with patch
end

function S_patch = construct_S_patch(f, theta_A, theta_B, h, n_xi, n_eta)
    l_theta = @(theta) [2*sin(theta/2), -sin(theta)];
    
    l_A = @(xi) (theta_B - theta_A)*xi + theta_A;
    nu = @(xi) [cos(l_A(xi)), cos(l_A(xi)/2)] ./ sqrt(cos(l_A(xi)).^2 + cos(l_A(xi)/2).^2);
    
    M_p_general = @(xi, eta, H) l_theta(l_A(xi)) + eta.*H.*nu(xi);%./nu_norm(xi);
    % H is a function of h and n_eta
    
    theta_diff = theta_B - theta_A;
    J = @(v, H) [theta_diff*cos(l_A(v(1))/2)-v(2)*H*theta_diff*sin(l_A(v(1))), H*cos(l_A(v(1))); -1*theta_diff*cos(l_A(v(1)))-v(2)*H*theta_diff/2*sin(l_A(v(1))/2), H*cos(l_A(v(1))/2)];
    
    xi_mesh = linspace(0, 1, n_xi+1);
    eta_mesh = linspace(0, 1, n_eta+1);
    
    [XI, ETA] = meshgrid(xi_mesh, eta_mesh);
    XY = M_p_general(XI(:), ETA(:), h*n_eta);
    
    f_XY = reshape(f(XY(:, 1), XY(:, 2)), [n_eta+1, n_xi+1]);
    
    S_patch = S_patch_obj(M_p_general, J, n_xi, n_eta, 0, 1, 0, 1, f_XY, h, nan);
end