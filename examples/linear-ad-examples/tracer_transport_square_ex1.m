clc; clear; close all;
%% Problem Parameters
% spatial discretization
h = 0.01;

% Peclet Number
Pe = 1000;

% initial condition and boundary conditions
f0 = @(x, y) 1/4*((x-1).^2 + y.^3 +1);

% pressure field (solution to Laplace equation)
p_field = @(x, y) (x.^2-y.^2);

t_final = 0.6;

%% FC Parameters
d = 7;
C = 27;
n_r = 1;

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

%% 2DFC of initial condition

x_mesh = (0:h:1)';
y_mesh = (0:h:1)';
n = length(x_mesh);

boundary_X = [x_mesh; ones(n, 1); flipud(x_mesh); zeros(n, 1)];
boundary_Y = [zeros(n, 1); y_mesh; ones(n, 1); flipud(y_mesh)];

[X, Y] = meshgrid(x_mesh, y_mesh);

R_c = R_cartesian_mesh_obj(-C*h, 1+C*h-h/2, -C*h, 1+C*h-h/2, h, boundary_X, boundary_Y);
R_c.in_interior = R_c.R_X > 0 & R_c.R_X < 1 & R_c.R_Y > 0 & R_c.R_Y < 1;
R_c.f_R = fcont_gram_blend_square(f0(X, Y), d, A, Q);

%% Pressure/flow computation
R_p = R_cartesian_mesh_obj(R_c.x_start, R_c.x_end-h/2, R_c.y_start, R_c.y_end-h/2, h, R_c.boundary_X, R_c.boundary_Y);
R_p.f_R = p_field(R_c.R_X, R_c.R_Y);

% Compute flow field with Spectral differentiation
u_x_exact = @(x, y) 2*x;
u_y_exact = @(x, y) -2*y;

u_x = u_x_exact(R_p.R_X, R_p.R_Y);
u_y = u_y_exact(R_p.R_X, R_p.R_Y);

%% Linear Advection-Diffusion
% time step size computed with Von-Neumann stability analysis
dt = 1/5*min(1.4*Pe/pi^2*h^2, 1.5/(pi*max(sqrt(u_x.^2+u_y.^2), [], 'all'))*h);

t_mesh = dt:dt:t_final;
size(t_mesh)

vid_time = 5
fr = 30;
cap_stepsize = floor(length(t_mesh)/fr/5)

figure(1);
v = VideoWriter('scatter_animation.mp4', 'MPEG-4'); % Create video file
v.FrameRate = fr; % Set frame rate
open(v);

num_figs = 4;
t_idx_figs = floor(linspace(1, length(t_mesh), num_figs));
curr_fig = 1;

for t_idx = 1:length(t_mesh)
    t = t_mesh(t_idx);
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
    
    R_c.f_R =  R_c.f_R + dt/6*(k1+2*k2+2*k3+k4);
    
    %imposing boundary conditions
    R_c.f_R(C+1:end-C, C+1) = f0(R_c.R_X(C+1:end-C, C+1), R_c.R_Y(C+1:end-C, C+1));
    R_c.f_R(C+1:end-C, end-C) = f0(R_c.R_X(C+1:end-C, end-C), R_c.R_Y(C+1:end-C, end-C));
    R_c.f_R(C+1, C+1:end-C) = f0(R_c.R_X(C+1, C+1:end-C), R_c.R_Y(C+1, C+1:end-C));
    R_c.f_R(end-C, C+1:end-C) = f0(R_c.R_X(end-C, C+1:end-C), R_c.R_Y(end-C, C+1:end-C));
    
%     if mod(t_idx, cap_stepsize) == 0
%         scatter(R_c.R_X(R_c.in_interior), R_c.R_Y(R_c.in_interior), 50,R_c.f_R(R_c.in_interior), 'filled'); % Color by z values
% %         caxis([0, 1])
%         colorbar; % Show color scale
%         xlabel('x'); ylabel('y'); title(['Concentration at t=', num2str(t)]);
%         drawnow;
%         
%         frame = getframe(gcf);
%         writeVideo(v, frame);
%     end
    if t_idx == t_idx_figs(curr_fig)
        subplot(1, num_figs, curr_fig)
         scatter(R_c.R_X(R_c.in_interior), R_c.R_Y(R_c.in_interior), 50,R_c.f_R(R_c.in_interior), 'filled'); % Color by z values
%         caxis([0, 1])
        colorbar; % Show color scale
        xlabel('x'); ylabel('y'); title(['Concentration at t=', num2str(t)]);
        curr_fig = curr_fig + 1;
    end
    
    
    R_c.f_R = fcont_gram_blend_square(R_c.f_R(C+1:end-C, C+1:end-C), d, A, Q);    
end
close(v);

set(figure(1), 'Position', [100, 100, 1024, 164])

function [fcont] = fcont_gram_blend_square(f_xy, d, A, Q)
    fcont1 = [fcont_gram_blend_S(f_xy, d, A, Q); f_xy(2:end, :)]; 
    fcont2 = [fcont1(1:end-1, :); flipud(fcont_gram_blend_S(flipud(fcont1), d, A, Q))];    
    fcont_lr = rot90(fcont2);
    fcont3 = [fcont_gram_blend_S(fcont_lr, d, A, Q); fcont_lr(2:end, :)]; 
    fcont4 = [fcont3(1:end-1, :); flipud(fcont_gram_blend_S(flipud(fcont3), d, A, Q))];
    
    fcont = rot90(fcont4, 3);
end