clc; clear; close all;

f = @(x, y) -cos(x) - sin(y);
u_boundary = @(x, y) cos(x)+sin(y);

d = 8;
C = 27;
n_r = 6;

M = d+3;

h = 0.00125;
h_tan = h;


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

figure;
hold on;
theta_mesh = linspace(0, 1, 100)';

l_1 = @(theta) theta;
l_2 = @(theta) zeros(size(theta));
l_1_prime = @(theta) ones(size(theta));
l_2_prime = @(theta) zeros(size(theta));
l_1_dprime = @(theta) zeros(size(theta));
l_2_dprime = @(theta) zeros(size(theta));

curve_seq.add_curve(l_1, l_2, l_1_prime, l_2_prime, l_1_dprime, l_2_dprime, 0, 1/3, 1/3, 2/3, 2/3, h_tan)

l_1 = @(theta) (1-theta)+theta*1/2;
l_2 = @(theta) theta*sqrt(3)/2;
l_1_prime = @(theta) -1/2*ones(size(theta));
l_2_prime = @(theta) sqrt(3)/2*ones(size(theta));
l_1_dprime = @(theta) zeros(size(theta));
l_2_dprime = @(theta) zeros(size(theta));

curve_seq.add_curve(l_1, l_2, l_1_prime, l_2_prime, l_1_dprime, l_2_dprime, 0, 1/3, 1/3, 2/3, 2/3, h_tan)

l_1 = @(theta) (1-theta)*1/2;
l_2 = @(theta) (1-theta)*sqrt(3)/2;
l_1_prime = @(theta) -1/2*ones(size(theta));
l_2_prime = @(theta) -sqrt(3)/2*ones(size(theta));
l_1_dprime = @(theta) zeros(size(theta));
l_2_dprime = @(theta) zeros(size(theta));

curve_seq.add_curve(l_1, l_2, l_1_prime, l_2_prime, l_1_dprime, l_2_dprime, 0, 1/3, 1/3, 2/3, 2/3, h_tan)

[u_num_mat, R] = poisson_solver(curve_seq, f, u_boundary, h, 1, 4, 1e-14, 1e-14, d, C, n_r, A, Q, M);
u_exact = u_boundary(R.R_X, R.R_Y);

max(abs(u_num_mat(R.in_interior)- u_exact(R.in_interior)), [], 'all')