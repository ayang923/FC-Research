clc; clear; close all;

f = @(x, y) -cos(x) - sin(y);
u_boundary = @(x, y) cos(x)+sin(y);

d = 4;
C = 27;
n_r = 6;

M = d+5;

h_tan = 0.004;
h = 0.002;

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

l_1 = @(theta) -theta+1;
l_2 = @(theta) -theta+1;
l_1_prime = @(theta) -1*ones(size(theta));
l_2_prime = @(theta) -1*ones(size(theta));
l_1_dprime = @(theta) zeros(size(theta));
l_2_dprime = @(theta) zeros(size(theta));

curve_seq.add_curve(l_1, l_2, l_1_prime, l_2_prime, l_1_dprime, l_2_dprime, 0, 1/3, 1/3, 2/3, 2/3, h_tan)

l_1 = @(theta) theta;
l_2 = @(theta) theta.^3;
l_1_prime = @(theta) ones(size(theta));
l_2_prime = @(theta) 3*theta.^2;
l_1_dprime = @(theta) zeros(size(theta));
l_2_dprime = @(theta) 6*theta;

curve_seq.add_curve(l_1, l_2, l_1_prime, l_2_prime, l_1_dprime, l_2_dprime, 0, 1/3, 1/3, 2/3, 2/3, h_tan)


[u_num_mat, R] = poisson_solver(curve_seq, f, u_boundary, h, 1, 4, 1e-14, 1e-14, d, C, n_r, A, Q, M);
u_exact = u_boundary(R.R_X, R.R_Y);

max(abs(u_num_mat(R.in_interior)- u_exact(R.in_interior)), [], 'all')