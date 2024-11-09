clc; clear; close all;

% 2DFC on boomerang domain
%% Setting Parameters
f = @(x, y) 4 + (1 + x.^2 + y.^2).*(sin(2.5*pi*x - 0.5) + cos(2*pi*y - 0.5));

alpha = 3/2;
beta = tan(alpha.*pi./2);
l_theta = @(theta) [-2./3.*sin(3./2.*theta), beta.*sin(theta)];

scale_factor = 8;

h_R = 0.005 / scale_factor;

n_C1 = round(0.2/(h_R))+1;
n_S_h = round(0.2/h_R);
h_S = h_R;

n_S_w = round((2*pi-0.5)/h_R);

eps_xi_eta = 1e-13;
eps_xy = 1e-13;

d = 7;

%% Constructing Boundary Patches
C1_patch = construct_C1_patch(f, 0.2, 2*pi-0.2, 0, -0.2, 2*pi+0.2, n_C1, n_C1, d, eps_xi_eta, eps_xy); % data associated with patch

window_patch_eta = construct_S_patch(f, 2*pi-0.05, 2*pi-0.25, h_S, n_S_h, d+5, eps_xi_eta, eps_xy);
window_patch_xi = construct_S_patch(f, 0.05, 0.25, h_S, n_S_h, d+5, eps_xi_eta, eps_xy);
S_patch = construct_S_patch(f, 0.25, 2*pi-0.25, h_S, n_S_w, d+5, eps_xi_eta, eps_xy);

save(['boomerang_data/patches_nC1bound', num2str(n_C1), '_d', num2str(d)], 'C1_patch', 'window_patch_xi', 'window_patch_eta', 'S_patch')

%% Computing Points of Interior Cartesian Mesh
R_x_bounds = [S_patch.Q_patch.x_min-h_R, S_patch.Q_patch.x_max+h_R];
R_y_bounds = [S_patch.Q_patch.y_min-h_R, S_patch.Q_patch.y_max+h_R];

% Computes boundary_XY values of domain 
boundary_XY = l_theta(transpose(linspace(0, 2*pi, 10000*scale_factor)));

% Constructs cartesian mesh object
R = R_cartesian_mesh_obj(R_x_bounds(1), R_x_bounds(2), R_y_bounds(1), R_y_bounds(2), h_R, boundary_XY(:, 1), boundary_XY(:, 2), true, nan);
R.fill_interior(f);

interior_idxs = R.interior_idxs;
f_R = R.f_R;

figure;
scatter3(R.R_X(R.interior_idxs), R.R_Y(R.interior_idxs), R.f_R(R.interior_idxs))

save(['boomerang_data/interior_mesh_nC1bound', num2str(n_C1)], 'R', 'h_R');
%% Problem specific patch construction
function S_patch = construct_S_patch(f, theta_A, theta_B, h, n_xi, n_eta, eps_xi_eta, eps_xy)
    alpha = 3/2;
    beta = tan(alpha.*pi./2);
    
    l_theta = @(theta) [-2./3.*sin(3./2.*theta), beta.*sin(theta)];
    nu_theta = @(theta) [beta.*cos(theta), cos(3./2.*theta)] ./ sqrt((beta.*cos(theta)).^2+(cos(3./2.*theta)).^2);
    l_A = @(xi) (theta_B - theta_A).*xi + theta_A;
    
    M_p_general = @(xi, eta, H) l_theta(l_A(xi)) + eta.*H.*nu_theta(l_A(xi));
    dM_pdxi = @(xi, eta, H) [
        -((-theta_A + theta_B) .* cos((3./2) .* (theta_A + (-theta_A + theta_B) .* xi))) ...
        - (beta .* eta .* H .* (-theta_A + theta_B) .* sin(theta_A + (-theta_A + theta_B) .* xi)) ...
        ./ sqrt(beta.^2 .* cos(theta_A + (-theta_A + theta_B) .* xi).^2 + cos((3./2) .* (theta_A + (-theta_A + theta_B) .* xi)).^2) ...
        - (beta .* eta .* H .* cos(theta_A + (-theta_A + theta_B) .* xi) ...
        .* (-2 .* beta.^2 .* (-theta_A + theta_B) .* cos(theta_A + (-theta_A + theta_B) .* xi) .* sin(theta_A + (-theta_A + theta_B) .* xi) ...
        - 3 .* (-theta_A + theta_B) .* cos((3./2) .* (theta_A + (-theta_A + theta_B) .* xi)) .* sin((3./2) .* (theta_A + (-theta_A + theta_B) .* xi)))) ...
        ./ (2 .* (beta.^2 .* cos(theta_A + (-theta_A + theta_B) .* xi).^2 + cos((3./2) .* (theta_A + (-theta_A + theta_B) .* xi)).^2).^(3./2)); ...
        
        beta .* (-theta_A + theta_B) .* cos(theta_A + (-theta_A + theta_B) .* xi) ...
        - (3 .* eta .* H .* (-theta_A + theta_B) .* sin((3./2) .* (theta_A + (-theta_A + theta_B) .* xi))) ...
        ./ (2 .* sqrt(beta.^2 .* cos(theta_A + (-theta_A + theta_B) .* xi).^2 + cos((3./2) .* (theta_A + (-theta_A + theta_B) .* xi)).^2)) ...
        - (eta .* H .* cos((3./2) .* (theta_A + (-theta_A + theta_B) .* xi)) ...
        .* (-2 .* beta.^2 .* (-theta_A + theta_B) .* cos(theta_A + (-theta_A + theta_B) .* xi) .* sin(theta_A + (-theta_A + theta_B) .* xi) ...
        - 3 .* (-theta_A + theta_B) .* cos((3./2) .* (theta_A + (-theta_A + theta_B) .* xi)) .* sin((3./2) .* (theta_A + (-theta_A + theta_B) .* xi)))) ...
        ./ (2 .* (beta.^2 .* cos(theta_A + (-theta_A + theta_B) .* xi).^2 + cos((3./2) .* (theta_A + (-theta_A + theta_B) .* xi)).^2).^(3./2));
    ];
    
    dM_pdeta = @(xi, eta, H) [
        (beta .* H .* cos(theta_A + (-theta_A + theta_B) .* xi)) ...
            ./ sqrt(beta.^2 .* cos(theta_A + (-theta_A + theta_B) .* xi).^2 + cos((3./2) .* (theta_A + (-theta_A + theta_B) .* xi)).^2); ...
        (H .* cos((3./2) .* (theta_A + (-theta_A + theta_B) .* xi))) ...
            ./ sqrt(beta.^2 .* cos(theta_A + (-theta_A + theta_B) .* xi).^2 + cos((3./2) .* (theta_A + (-theta_A + theta_B) .* xi)).^2);
    ];
    
    J = @(v, H) [dM_pdxi(v(1), v(2), H), dM_pdeta(v(1), v(2), H)];
    
    S_patch = S_patch_obj(M_p_general, J, h, eps_xi_eta, eps_xy, n_xi, n_eta, nan);

    [X, Y] = S_patch.Q_patch.xy_mesh;
    S_patch.Q_patch.f_XY = f(X, Y);
end

function C1_patch = construct_C1_patch(f, theta_A, theta_B, theta_C, theta_D, theta_E, n_xi, n_eta, d, eps_xi_eta, eps_xy)
    alpha = 3/2;
    beta = tan(alpha.*pi./2);
    
    l_theta = @(theta) [-2./3.*sin(3./2.*theta), beta.*sin(theta)];
    l_A = @(xi) (theta_D - theta_A)*xi + theta_A;
    l_B = @(eta) (theta_E - theta_B)*eta + theta_B;
    
    M_p = @(xi, eta) l_theta(l_A(xi)) + l_theta(l_B(eta)) - l_theta(theta_C);
    
    dM_pdxi = @(xi, eta) [-((-theta_A+theta_D)*cos(3/2*(theta_A+(-theta_A+theta_D)*xi))); beta*(-theta_A+theta_D)*cos(theta_A+(-theta_A+theta_D)*xi)];
    dM_pdeta = @(xi, eta) [-((-theta_B+theta_E)*cos(3/2*(theta_B+(-theta_B+theta_E)*eta))); beta*(-theta_B+theta_E)*cos(theta_B+(-theta_B+theta_E)*eta)];

    J = @(v) [dM_pdxi(v(1), v(2)), dM_pdeta(v(1), v(2))];
    C1_patch = C1_patch_obj(M_p, J, eps_xi_eta, eps_xy, n_xi, n_eta, d, nan, nan);
    
    [X_L, Y_L] = C1_patch.L.xy_mesh;
    C1_patch.L.f_XY = f(X_L, Y_L);
    
    [X_W, Y_W] = C1_patch.W.xy_mesh;
    C1_patch.W.f_XY = f(X_W, Y_W);
end