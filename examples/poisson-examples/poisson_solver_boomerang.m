clc; clear; close all;

f = @(x, y) -cos(x) - sin(y);
u_boundary = @(x, y) cos(x)+sin(y);

d = 4;
C = 27;
n_r = 6;

M = d+3;
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

alpha = 3/2;
beta = tan(alpha.*pi./2);

l_1 = @(theta) -2/3*sin(theta*3*pi);
l_2 = @(theta) beta*sin(theta*2*pi);
l_1_prime = @(theta) -2*pi*cos(theta*3*pi);
l_2_prime = @(theta) 2*pi*beta*cos(theta*2*pi);
l_1_dprime = @(theta) 6*pi^2*sin(theta*3*pi);
l_2_dprime = @(theta) -4*pi^2*beta*sin(theta*2*pi);

h = 0.005;

curve_seq = Curve_seq_obj();
curve_seq.add_curve(l_1, l_2, l_1_prime, l_2_prime, l_1_dprime, l_2_dprime, 0, 1/10, 1/10, 0, 0, h);

[u_num_mat, R] = poisson_solver(curve_seq, f, u_boundary, h, 1, 4, 1e-14, 1e-14, d, C, n_r, A, Q, M);

u_exact = u_boundary(R.R_X, R.R_Y);

max(abs(u_num_mat(R.in_interior)- u_exact(R.in_interior)), [], 'all')