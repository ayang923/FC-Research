clc; clear; close all;

%% Problem Parameters
% spatial discretization
h = 0.01;

% Peclet Number
Pe = 100;

% initial condition and boundary conditions
f = @(x, y) 1/4*((x-1).^2 + y.^3 +1);

% pressure field (solution to Laplace equation)
p_field = @(x, y) (x.^2-y.^2)./4;

t_final = 1;

%% Domain Parametrization
% setting up domain
l_1 = @(theta) 2*sin(theta*pi);
l_2 = @(theta) -sin(theta*2*pi);
l_1_prime = @(theta) 2*pi*cos(theta*pi);
l_2_prime = @(theta) -2*pi*cos(theta*2*pi);
l_1_dprime = @(theta) -2*pi^2*sin(theta*pi);
l_2_dprime = @(theta) 4*pi^2*sin(theta*2*pi);

%% FC Parameters
d = 4;
C = 27;
n_r = 6;
M = d+3;

n_frac_C = 1/10;

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

curve_seq = Curve_seq_obj();
curve_seq.add_curve(l_1, l_2, l_1_prime, l_2_prime, l_1_dprime, l_2_dprime, 0, n_frac_C, n_frac_C, 0, 0, h);

eps_xi_eta = 1e-13;
eps_xy = 1e-13;

%% 2DFC of initial condition, also generates Cartesian mesh and patch meshes
[R_c, interior_patches, fc_patches] = FC2D(f, h, curve_seq, eps_xi_eta, eps_xy, d, C, n_r, A, Q, C, A, Q, M);

figure;
surf(R_c.R_X, R_c.R_Y, R_c.f_R)

%wrapping R_c as a q_patch to use locally compute
R_c_q_patch = Q_patch_obj(@(xi, eta) [xi, eta], nan, eps_xy, eps_xy, R_c.n_x, R_c.n_y, R_c.x_start, R_c.x_end, R_c.y_start, R_c.y_end, R_c.f_R);
%saving meshes and partition of unity
interior_patch_meshes = curve_seq.construct_patches(@(x, y) x*0+1, d, eps_xi_eta, eps_xy);
%% Pressure/flow computation

% Pressure mesh for 2DFC
R_p = R_cartesian_mesh_obj(R_c.x_start, R_c.x_end-h/2, R_c.y_start, R_c.y_end-h/2, h, R_c.boundary_X, R_c.boundary_Y);
p_interior = p_field(R_c.R_X(R_c.in_interior), R_c.R_Y(R_c.in_interior));

% Constructing patches for pressure
p_patches = cell(length(interior_patch_meshes));
for i = 1:curve_seq.n_curves
    S_patch = interior_patch_meshes{2*i-1}.copy_patch;
    C_patch = interior_patch_meshes{2*i}.copy_patch;
    
    [X_S, Y_S] = S_patch.Q.xy_mesh;
    [X_L, Y_L] = C_patch.L.xy_mesh;
    [X_W, Y_W] = C_patch.W.xy_mesh;
    
    S_patch.Q.f_XY = p_field(X_S, Y_S).*interior_patch_meshes{2*i-1}.Q.f_XY;
    C_patch.L.f_XY = p_field(X_L, Y_L).*interior_patch_meshes{2*i}.L.f_XY;
    C_patch.W.f_XY = p_field(X_W, Y_W).*interior_patch_meshes{2*i}.W.f_XY;
    
    p_patches{2*i-1} = S_patch;
    p_patches{2*i} = C_patch;
end

% figure;
% surf(X_S, Y_S, S_patch.Q.f_XY)
% hold on;
% surf(X_W, Y_W, C_patch.W.f_XY)
% surf(X_L, Y_L, C_patch.L.f_XY)
% drawnow;
% hold off;


% Performs 2DFC on pressure mesh
FC2D_patches(p_patches, R_p, p_interior, curve_seq, d, C, n_r, A, Q, C, A, Q, M);

% Compute flow field with Spectral differentiation
u_x_exact = @(x, y) x/2;
u_y_exact = @(x, y) -y/2;

[u_x, u_y] = R_p.grad(R_p.f_R);

u_div = R_p.div(u_x, u_y);
u_lap = R_p.lap(R_p.f_R);
max(abs(u_div(R_p.in_interior)))
max(abs(u_lap(R_p.in_interior)))
% u_x = R_u_x.f_R;
% u_y = R_u_y.f_R;

u_x = u_x_exact(R_p.R_X, R_p.R_Y);
u_y = u_x_exact(R_p.R_X, R_p.R_Y);

%% Linear Advection-Diffusion
% time step size computed with Von-Neumann stability analysis
dt = 1/5*min(1.4*Pe/pi^2*h^2, 1.5/(pi*max(sqrt(u_x.^2+u_y.^2), [], 'all'))*h)

t_mesh = dt:dt:t_final;
size(t_mesh)

hold off;
figure;
% scatter(R_c.R_X(R_c.in_interior), R_c.R_Y(R_c.in_interior), 50, R_c.f_R(R_c.in_interior), 'filled'); % Color by z values
% caxis([0, 1]);
% colorbar; % Show color scale
surf(R_c.R_X, R_c.R_Y, R_c.f_R)
xlabel('x'); ylabel('y'); title('Concentration vs Time');

drawnow;

v = VideoWriter('scatter_animation.mp4', 'MPEG-4'); % Create video file
v.FrameRate = 30; % Set frame rate
open(v);

for t = t_mesh
    tic;
    y_k1 = R_c.f_R;
    [lap_c1] = R_c.lap(y_k1);
    [div_uc1] = R_c.div(u_x.*y_k1, u_y.*y_k1);
    k1 = -div_uc1 + 1/Pe*lap_c1;

    y_k2 = R_c.f_R+dt*k1/2;
    [lap_c2] = R_c.lap(y_k2);
    [div_uc2] = R_c.div(u_x.*y_k2, u_y.*y_k2);
    k2 = -div_uc2 + 1/Pe*lap_c2;

    y_k3 = R_c.f_R+dt*k2/2;
    [lap_c3] = R_c.lap(y_k3);
    [div_uc3] = R_c.div(u_x.*y_k3, u_y.*y_k3);
    k3 = -div_uc3 + 1/Pe*lap_c3;

    y_k4 = R_c.f_R+dt*k3;
    [lap_c4] = R_c.lap(y_k4);
    [div_uc4] = R_c.div(u_x.*y_k4, u_y.*y_k4);
    k4 = -div_uc4 + 1/Pe*lap_c4;
    
    f_np1 =  R_c.f_R + dt/6*(k1+2*k2+2*k3+k4);
%     f_np1 =  R_c.f_R;
    R_c_q_patch.f_XY = f_np1;
        
    for i = 1:curve_seq.n_curves
        S_patch = interior_patches{2*i-1};
        C_patch = interior_patches{2*i}.copy_patch;

        [X_S, Y_S] = S_patch.Q.xy_mesh;
        [X_L, Y_L] = C_patch.L.xy_mesh;
        [X_W, Y_W] = C_patch.W.xy_mesh;

        f_XY_S = zeros(size(X_S));
%         f_XY_S(1, :) = f(X_S(1, :), Y_S(1, :)); % impose boundary condition
        for patch_i = 1:size(f_XY_S, 1)
            for patch_j = 1:size(f_XY_S, 2)
                f_XY_S(patch_i, patch_j) = R_c_q_patch.locally_compute(X_S(patch_i, patch_j), Y_S(patch_i, patch_j), M);
            end
        end

        f_XY_L = zeros(size(X_L));
        f_XY_W = zeros(size(X_W));
        if isa(C_patch, 'C2_patch_obj')
%             f_XY_L(:, 1) = f(X_L(:, 1), Y_L(:, 1));
            for patch_i = 1:size(f_XY_L, 1)
                for patch_j = 1:size(f_XY_L, 2)
                    f_XY_L(patch_i, patch_j) = R_c_q_patch.locally_compute(X_L(patch_i, patch_j), Y_L(patch_i, patch_j), M);
                end
            end

%             f_XY_W(1, :) = f(X_W(1, :), Y_W(1, :));
            for patch_i = 1:size(f_XY_W, 1)
                for patch_j = 1:size(f_XY_W, 2)
                    f_XY_W(patch_i, patch_j) = R_c_q_patch.locally_compute(X_W(patch_i, patch_j), Y_W(patch_i, patch_j), M);
                end
            end
        else
            f_XY_L(1, :) = f(X_L(1, :), Y_L(1, :));
            for patch_i = 2:size(f_XY_L, 1)
                for patch_j = 1:size(f_XY_L, 2)
                    f_XY_L(patch_i, patch_j) = R_c_q_patch.locally_compute(X_L(patch_i, patch_j), Y_L(patch_i, patch_j), M);
                end
            end

            f_XY_W(:, 1) = f(X_W(:, 1), Y_W(:, 1));
            for patch_i = 1:size(f_XY_W, 1)
                for patch_j = 2:size(f_XY_W, 2)
                    f_XY_W(patch_i, patch_j) = R_c_q_patch.locally_compute(X_W(patch_i, patch_j), Y_W(patch_i, patch_j), M);
                end
            end
        end
        S_patch.Q.f_XY = f_XY_S.*interior_patch_meshes{2*i-1}.Q.f_XY;
        C_patch.L.f_XY = f_XY_L.*interior_patch_meshes{2*i}.L.f_XY;
        C_patch.W.f_XY = f_XY_W.*interior_patch_meshes{2*i}.W.f_XY;
    end
    
    
%     scatter(R_c.R_X, R_c.R_Y(R_c.in_interior), 50, f_np1(R_c.in_interior), 'filled'); % Color by z values
%     caxis([0, 1])
%     colorbar; % Show color scale
    xlabel('x'); ylabel('y'); title('Concentration vs Time');
%     drawnow;
    
    frame = getframe(gcf);
    writeVideo(v, frame);

    R_c.f_R = 0*R_c.f_R;
    FC2D_patches(interior_patches, R_c, f_np1(R_c.in_interior), curve_seq, d, C, n_r, A, Q, C, A, Q, M);
    disp(num2str(t))
    
    surf(R_c.R_X, R_c.R_Y, R_c.f_R)
    drawnow;
    toc;
end

close(v)