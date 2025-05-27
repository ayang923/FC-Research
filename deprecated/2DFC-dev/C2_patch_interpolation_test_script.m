% clc; clear; close all
% rng("default")
% 
% f = @(x, y) 4 + (1 + x.^2 + y.^2).*(sin(2.5*pi*x - 0.5) + cos(2*pi*y - 0.5));
% l_theta = @(theta) [2*sin(theta/2), -sin(theta)];
% 
% n_lst = 20 * 2.^(5:6);
% h_lst = 0.029 * (1/2).^(5:6);
% 
% C2_patches = cell(length(n_lst), 1);
% C2_norms = cell(length(h_lst), 1);
% 
% % Construct FC Patches
% for i=1:length(n_lst)
%     n_C2 = n_lst(i);
%     h_S = h_lst(i);
%     n_S = n_C2; % 10;%ceil(21/20*n_C2);
%     
%     C2_patch = construct_C2_patch(f, 0.4, 2*pi-0.4, 0, n_C2, n_C2);
%     window_patch_xi = construct_S_patch(f, 0.1, 0.6, h_S, n_S, n_S);
%     window_patch_eta = construct_S_patch(f, 2*pi-0.1, 2*pi-0.6, h_S, n_S, n_S);
%     
%     [C2_norm, xi_norm, eta_norm] = C2_patch.compute_phi_normalization(window_patch_xi, window_patch_eta);
%     
%     C2_patches{i} = C2_patch;
%     C2_norms{i} = C2_norm;
%     
% end
% 
% save("C2_patches.mat", "C2_patches", "C2_norms")
% 
% disp("patches constructed")

clc; clear; close all;



n_lst = 20 * 2.^(5:6);
h_lst = 0.029 * (1/2).^(5:6);

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

finest_C2_patch.f_XY = f(XI_err, ETA_err); %finest_C2_patch.phi(XI_err, ETA_err); % ./ C2_norms{length(n_lst)};
finest_C2_fcont_patch = finest_C2_patch.FC(C, d, A, Q, nan);
test_patch_i = 1;
% FC and FFT
C2_patch = C2_patches{test_patch_i};
C2_norm = C2_norms{test_patch_i};

[XI, ETA] = C2_patch.xi_eta_mesh();
% C2_patch.f_XY = C2_patch.f_XY .* C2_patch.phi(XI, ETA);

C2_patch.f_XY = f(XI, ETA); %C2_patch.phi(XI, ETA); %./C2_norm;
C2_fcont_patch = C2_patch.FC(C, d, A, Q, nan);

[XI_fcont, ETA_fcont] = C2_fcont_patch.xi_eta_mesh();

[bound_C2_X, bound_C2_Y] = C2_patch.boundary_mesh_xy();

R_x_bounds = [C2_fcont_patch.x_min, C2_fcont_patch.x_max];
R_y_bounds = [C2_fcont_patch.y_min, C2_fcont_patch.y_max];

disp("creating meshes")
R = R_cartesian_mesh_obj(R_x_bounds(1), R_x_bounds(2), R_y_bounds(1), R_y_bounds(2), h_lst(1), bound_C2_X, bound_C2_Y);
R_err = R_cartesian_mesh_obj(R_x_bounds(1), R_x_bounds(2), R_y_bounds(1), R_y_bounds(2), h_lst(1)./2, bound_C2_X, bound_C2_Y);

disp("meshes created");

R.interpolate_patch(C2_fcont_patch, d+3, false);
R_err.interpolate_patch(C2_fcont_patch, d+3, false);
disp("interpolation done")

save("R_err.mat", "R_err")
save("R.mat", "R");

clc; clear;

load("R.mat");
load("R_err.mat");

disp("start FFT")
R.compute_fc_coeffs();
[R_X_err, R_Y_err, f_interpolation, interior_idx] = R.ifft_interpolation(R.h./2);
disp("FFT finished")

err = abs(R_err.f_R(R_err.interior_idxs)-f_interpolation(interior_idx));
max(err)

figure;
s = surf(R_err.R_X, R_err.R_Y, R_err.f_R);
s.EdgeColor = 'none';

err_mesh = R.y_start:R.h/2:R.y_end+R.h/2;
N = length(err_mesh);

load("C2_patches")

errs = zeros(size(R.f_R, 2), 1);
for i = 1:size(R.f_R, 2)
    fc_coeffs = fftshift(1/R.n_y * fft(R.f_R(:, i)));

    padded_fc_coeffs = [zeros(floor((N-R.n_y)/2), 1); fc_coeffs; zeros(ceil((N-R.n_y)/2), 1)];
    f_numeric = N*real(ifft(ifftshift(padded_fc_coeffs)));
    
    if 2*i-1 < size(R_err.f_R, 2)
        errs(i) = max(abs(f_numeric(1:end-1) - R_err.f_R(:, 2*i-1)));
    end
end

figure;
plot(errs)

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
    
    M_p_general = @(xi, eta, H) l_theta(l_A(xi)) + eta.*H.*nu(xi);
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

function f_y = f(x, y)
    r = sqrt(x.^2 + y.^2);
    f_y = zeros(size(x));
    f_y(r < 0.56) = exp(1./(1/0.56*r(r < 0.56).^2-0.56));
    f_y(r >= 0.56) = 0;
end