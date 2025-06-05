clc; clear; close all;
load('geo_data_poly.mat')

p = 2;
f = @(x, y) x.^2 - y.^2;

v = @(s) (1/p-1/2)*(1-2*s).^3+1/p*(2*s-1)+1/2;
v_prime = @(s) 2/p-6*(1/p-1/2)*(1-2*s).^2;

w = @(s) v(s).^p./(v(s).^p+v(1-s).^p);
w_prime = @(s) p*((v_prime(s).*v(s).^(p-1))./(v(s).^p+v(1-s).^p)-(v(s).^(p-1).*v_prime(s)-v(1-s).^(p-1).*v_prime(1-s)).*v(s).^p./(v(s).^p+v(1-s).^p).^2);

[A, b, n_total, curve_n, start_idx, end_idx] = construct_A_b(1, f, curve_seq, w, w_prime);
phi_j = A\b;

gr_phi = zeros(n_total, 1);
curr = curve_seq.first_curve;
for i = 1:curve_seq.n_curves
    ds = 1/curve_n(i);
    s_mesh = linspace(0, 1-ds, curve_n(i))';

    gr_phi(start_idx(i):end_idx(i)) = w_prime(s_mesh).*phi_j(start_idx(i):end_idx(i));
    curr = curr.next_curve;
end


u_num = @(x, y) u_num_global(x, y, gr_phi, curve_seq, start_idx, end_idx, curve_n, w);

% M = 7;

R = 2;
fft_coeffs = fftshift(fft(gr_phi))/n_total;
padded_fft_coeffs = [zeros(floor((n_total*R-n_total)/2), 1); fft_coeffs; zeros(ceil((n_total*R-n_total)/2), 1)];
gr_phi_R = R*n_total*real(ifft(ifftshift(padded_fft_coeffs)));

[A_R, b_R] = construct_A_b(R, f, curve_seq, w, w_prime);
u_num_b(1, 10, gr_phi, curve_seq, start_idx, end_idx, curve_n, w)


function [A, b, n_total, curve_n, start_idx, end_idx] = construct_A_b(R, f, curve_seq, w, w_prime)
    curve_n = zeros(curve_seq.n_curves, 1);
    curr = curve_seq.first_curve;
    for i = 1:curve_seq.n_curves
        curve_n(i) = (curr.n-1)*R;
        curr = curr.next_curve;
    end

    n_total = sum(curve_n);
    start_idx = cumsum([1, curve_n(1:end-1)'])';
    end_idx = start_idx+curve_n-1;

    K_boundary = @(theta_1, theta_2, curve_1, curve_2) ...
        -1 / (2 * pi) * ( ...
            (curve_1.l_1(theta_1) - curve_2.l_1(theta_2)) .* curve_2.l_2_prime(theta_2) - ...
            (curve_1.l_2(theta_1) - curve_2.l_2(theta_2)) .* curve_2.l_1_prime(theta_2) ...
        ) ./ ( ...
            sqrt(curve_2.l_1_prime(theta_2).^2 + curve_2.l_2_prime(theta_2).^2) .* ...
            ( ...
                (curve_1.l_1(theta_1) - curve_2.l_1(theta_2)).^2 + ...
                (curve_1.l_2(theta_1) - curve_2.l_2(theta_2)).^2 ...
            ) ...
        );

    K_boundary_same_point = @(theta, curve) ...
        -1 / (4 * pi) * ( ...
            curve.l_2_prime(theta) .* curve.l_1_dprime(theta) - ...
            curve.l_1_prime(theta) .* curve.l_2_dprime(theta) ...
        ) ./ ( ...
            (curve.l_1_prime(theta).^2 + curve.l_2_prime(theta).^2).^(3/2) ...
        );

    A = zeros(n_total, n_total);
    b = zeros(n_total, 1);


    curr_y = curve_seq.first_curve;
    for i = 1:curve_seq.n_curves
        ds_y = 1/curve_n(i);
        s_mesh_y = linspace(0, 1-ds_y, curve_n(i))';
        theta_mesh_y = w(s_mesh_y);
        b(start_idx(i):end_idx(i)) = f(curr_y.l_1(theta_mesh_y), curr_y.l_2(theta_mesh_y));

        curr_x = curve_seq.first_curve;
        for j = 1:curve_seq.n_curves
            ds_x = 1/curve_n(j);
            s_mesh_x = linspace(0, 1-ds_x, curve_n(j))';
            theta_mesh_x = w(s_mesh_x);
            [theta_2, theta_1] = meshgrid(theta_mesh_y, theta_mesh_x);

            A_local = w_prime(s_mesh_y').*K_boundary(theta_1, theta_2, curr_x, curr_y).*sqrt(curr_y.l_1_prime(theta_mesh_y').^2+curr_y.l_2_prime(theta_mesh_y').^2)*ds_y;

            if i == j
                msk_diagonals = logical(diag(ones(curve_n(j), 1)));
                A_local(msk_diagonals) = w_prime(s_mesh_y).*K_boundary_same_point(theta_mesh_y, curr_x).*sqrt(curr_x.l_1_prime(theta_mesh_x).^2+curr_x.l_2_prime(theta_mesh_x).^2)*ds_x + 1/2;
            end

            A(start_idx(j):end_idx(j), start_idx(i):end_idx(i)) = A_local;

            curr_x = curr_x.next_curve;
        end

        curr_y = curr_y.next_curve;
    end
end

function u_num_b = u_num_b(curve_idx, s_idx, gr_phi, curve_seq, start_idx, end_idx, curve_n, w)
     K_boundary = @(theta_1, theta_2, curve_1, curve_2) ...
        -1 / (2 * pi) * ( ...
            (curve_1.l_1(theta_1) - curve_2.l_1(theta_2)) .* curve_2.l_2_prime(theta_2) - ...
            (curve_1.l_2(theta_1) - curve_2.l_2(theta_2)) .* curve_2.l_1_prime(theta_2) ...
        ) ./ ( ...
            sqrt(curve_2.l_1_prime(theta_2).^2 + curve_2.l_2_prime(theta_2).^2) .* ...
            ( ...
                (curve_1.l_1(theta_1) - curve_2.l_1(theta_2)).^2 + ...
                (curve_1.l_2(theta_1) - curve_2.l_2(theta_2)).^2 ...
            ) ...
        );

    K_boundary_same_point = @(theta, curve) ...
        -1 / (4 * pi) * ( ...
            curve.l_2_prime(theta) .* curve.l_1_dprime(theta) - ...
            curve.l_1_prime(theta) .* curve.l_2_dprime(theta) ...
        ) ./ ( ...
            (curve.l_1_prime(theta).^2 + curve.l_2_prime(theta).^2).^(3/2) ...
        );
    
    compute_curve = nan;
    curr = curve_seq.first_curve;
    for i = 1:size(start_idx)
        if i == curve_idx
            compute_curve = curr;
            break;
        end
        curr = curr.next_curve;
    end
    
    curr = curve_seq.first_curve;
    u_num_b = 0;
    for i = 1:size(start_idx)
        ds = 1/curve_n(i);
        s_mesh = linspace(0, 1-ds, curve_n(i))';
        theta_mesh = w(s_mesh);

        int_vec = transpose(K_boundary(w(s_mesh(s_idx)), theta_mesh, compute_curve, curr).*gr_phi(start_idx(i):end_idx(i)).*sqrt(curr.l_1_prime(theta_mesh).^2+curr.l_2_prime(theta_mesh).^2))* ds;
        if curve_idx == i
            int_vec(s_idx) = K_boundary_same_point(w(s_mesh(s_idx)), curr);
        end

        u_num_b = u_num_b + int_vec * ones(curve_n(i), 1);
        curr = curr.next_curve;
    end
end


function u_num = u_num_global(x, y, gr_phi_j, curve_seq, start_idx, end_idx, curve_n, w)
    K_general = @(x, theta, curve) ...
    -1 / (2 * pi) * ( ...
        (x(1) - curve.l_1(theta)) .* curve.l_2_prime(theta) - ...
        (x(2) - curve.l_2(theta)) .* curve.l_1_prime(theta) ...
    ) ./ ( ...
        sqrt(curve.l_1_prime(theta).^2 + curve.l_2_prime(theta).^2) .* ...
        ( ...
            (x(1) - curve.l_1(theta)).^2 + ...
            (x(2) - curve.l_2(theta)).^2 ...
        ) ...
    );

   curr = curve_seq.first_curve;
   u_num = 0;
   for i = 1:size(start_idx)
        ds = 1/curve_n(i);
        s_mesh = linspace(0, 1-ds, curve_n(i))';
        theta_mesh = w(s_mesh);

        u_num =  u_num + transpose(K_general([x; y], theta_mesh, curr).*gr_phi_j(start_idx(i):end_idx(i)).*sqrt(curr.l_1_prime(theta_mesh).^2+curr.l_2_prime(theta_mesh).^2)) * ones(curve_n(i), 1) * ds;
        curr = curr.next_curve;
    end
end