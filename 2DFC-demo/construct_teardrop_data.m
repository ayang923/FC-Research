clc; clear; close all;

%% Setting parameters
f = @(x, y) 4 + (1 + x.^2 + y.^2).*(sin(2.5*pi*x - 0.5) + cos(2*pi*y - 0.5));
l_theta = @(theta) [2*sin(theta/2), -sin(theta)];

scale_factor = 2;

h_R = 0.01 / scale_factor;

n_C2 = round(0.8/(h_R*2));
n_S_h = round(0.4/h_R);
h_S = h_R;

n_S_w = round((2*pi-1)/h_R);

eps_xi_eta = 1e-13;
eps_xy = 1e-13;

d = 7;

%% Constructing Boundary Patches
C2_patch = construct_C2_patch(f, 0.4, 2*pi-0.4, 0, n_C2, n_C2, eps_xi_eta, eps_xy); % data associated with patch

window_patch_xi = construct_S_patch(f, 0.1, 0.5, h_S, n_S_h, d+4, eps_xi_eta, eps_xy);
window_patch_eta = construct_S_patch(f, 2*pi-0.1, 2*pi-0.5, h_S, n_S_h, d+4, eps_xi_eta, eps_xy);
S_patch = construct_S_patch(f, 0.5, 2*pi-0.5, h_S, n_S_w, d+4, eps_xi_eta, eps_xy);

save(['teardrop_data/patches_nC2', num2str(n_C2), '_d', num2str(d)], 'C2_patch', 'window_patch_xi', 'window_patch_eta', 'S_patch')

%% Computing Points of Interior Cartesian Mesh
R_x_bounds = [C2_patch.x_min-h_R, S_patch.x_max+h_R];
R_y_bounds = [S_patch.y_min-h_R, S_patch.y_max+h_R];

% Computes boundary_XY values of domain 
boundary_XY = l_theta(transpose(linspace(0, 2*pi, 10000*scale_factor)));

% Constructs cartesian mesh object
R = R_cartesian_mesh_obj(R_x_bounds(1), R_x_bounds(2), R_y_bounds(1), R_y_bounds(2), h_R, boundary_XY(:, 1), boundary_XY(:, 2), true, nan);
R.fill_interior(f);

interior_idxs = R.interior_idxs;
f_R = R.f_R;

figure;
scatter3(R.R_X(R.interior_idxs), R.R_Y(R.interior_idxs), R.f_R(R.interior_idxs))

save(['teardrop_data/interior_mesh_nC2', num2str(n_C2)], 'R', 'h_R');

%% Problem specific patch construction
function C2_patch = construct_C2_patch(f, theta_A, theta_B, theta_C, n_xi, n_eta, eps_xi_eta, eps_xy)
    l_theta = @(theta) [2*sin(theta/2), -sin(theta)];
    
    l_A = @(xi) xi*(theta_A - theta_C) + theta_C;
    l_B = @(eta) eta*(theta_B - theta_C - 2*pi) + theta_C + 2*pi; % account for wrap around
    
    M_p = @(xi, eta) l_theta(l_A(xi)) + l_theta(l_B(eta)) - l_theta(theta_C);
    J = @(v) [theta_A*cos(v(1)*theta_A/2) (theta_B-2*pi)*cos((v(2)*(theta_B-2*pi)+2*pi)/2); -theta_A*cos(v(1)*theta_A) -(theta_B-2*pi)*cos(v(2)*(theta_B-2*pi)+2*pi)];
    
    xi_mesh = linspace(0, 1, n_xi);
    eta_mesh = linspace(0, 1, n_eta);
    
    [XI, ETA] = meshgrid(xi_mesh, eta_mesh);
    XY = M_p(XI(:), ETA(:));
    
    f_XY = reshape(f(XY(:, 1), XY(:, 2)), [n_eta, n_xi]);
    
    C2_patch = C2_patch_obj(M_p, J, eps_xi_eta, eps_xy, n_xi, n_eta, 0, 1, 0, 1, f_XY, nan); % data associated with patch
end

function S_patch = construct_S_patch(f, theta_A, theta_B, h, n_xi, n_eta, eps_xi_eta, eps_xy)
    l_theta = @(theta) [2*sin(theta/2), -sin(theta)];
    
    l_A = @(xi) (theta_B - theta_A)*xi + theta_A;
    nu = @(xi) [cos(l_A(xi)), cos(l_A(xi)/2)] ./ sqrt(cos(l_A(xi)).^2 + cos(l_A(xi)/2).^2);
    
    M_p_general = @(xi, eta, H) l_theta(l_A(xi)) + eta.*H.*nu(xi);%./nu_norm(xi);
    % H is a function of h and n_eta
    
    theta_diff = theta_B - theta_A;
    J = @(v, H) [theta_diff*cos(l_A(v(1))/2)-v(2)*H*theta_diff*sin(l_A(v(1))), H*cos(l_A(v(1))); -1*theta_diff*cos(l_A(v(1)))-v(2)*H*theta_diff/2*sin(l_A(v(1))/2), H*cos(l_A(v(1))/2)];
    
    xi_mesh = linspace(0, 1, n_xi);
    eta_mesh = linspace(0, 1, n_eta);
    
    [XI, ETA] = meshgrid(xi_mesh, eta_mesh);
    XY = M_p_general(XI(:), ETA(:), h*(n_eta-1));
    
    f_XY = reshape(f(XY(:, 1), XY(:, 2)), [n_eta, n_xi]);
    
    S_patch = S_patch_obj(M_p_general, J, eps_xi_eta, eps_xy, n_xi, n_eta, 0, 1, 0, 1, f_XY, h, nan);
end