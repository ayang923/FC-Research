clc; clear; close all;

f = @(x, y) 4 + (1 + x.^2 + y.^2).*(sin(2.5*pi*x - 0.5) + cos(2*pi*y - 0.5));

d = 4;
C_S = 27;
C_C = 27;
n_r = 6;

M = d+5;

h = 0.005;
h_tan = h;


if(exist(['FC_data/A_d',num2str(d),'_C', num2str(C_S), '_r', num2str(n_r), '.mat']) == 0 | ...
   exist(['FC_data/Q_d',num2str(d),'_C', num2str(C_S),  '_r', num2str(n_r), '.mat']) == 0)
    disp('FC data not found. Generating FC operators... \n');
    generate_bdry_continuations(d, C_S, C_S, 12, 20, 4, ...
        256, n_r);
end

load(['FC_data/A_d',num2str(d),'_C', num2str(C_S), '_r', num2str(n_r), '.mat']);
load(['FC_data/Q_d',num2str(d),'_C', num2str(C_S),  '_r', num2str(n_r), '.mat']);

A_S = double(A);
Q_S = double(Q);

if(exist(['FC_data/A_d',num2str(d),'_C', num2str(C_C), '_r', num2str(n_r), '.mat']) == 0 | ...
   exist(['FC_data/Q_d',num2str(d),'_C', num2str(C_C),  '_r', num2str(n_r), '.mat']) == 0)
    disp('FC data not found. Generating FC operators... \n');
    generate_bdry_continuations(d, C_C, C_C, 12, 20, 4, ...
        256, n_r);
end

load(['FC_data/A_d',num2str(d),'_C', num2str(C_C), '_r', num2str(n_r), '.mat']);
load(['FC_data/Q_d',num2str(d),'_C', num2str(C_C),  '_r', num2str(n_r), '.mat']);

A_C = double(A);
Q_C = double(Q);
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

[R, fc_patches] = FC2D(f, h, curve_seq, 5e-13, 1e-13, d, C_S, n_r, A_S, Q_S, C_C, A_C, Q_C, M);