function [err_] = laplace_solver_teardrop(f_boundary, x_start, x_end, y_start, y_end, h, interior_patches,  err_tol)
%      % teardrop
    l_1 = @(theta) 2*sin(theta*pi);
    l_2 = @(theta) -sin(theta*2*pi);
    l_1_prime = @(theta) 2*pi*cos(theta*pi);
    l_2_prime = @(theta) -2*pi*cos(theta*2*pi);
    l_1_dprime = @(theta) -2*pi^2*sin(theta*pi);
    l_2_dprime = @(theta) 4*pi^2*sin(theta*2*pi);

%     l_1 = @(theta) cos(2*pi*theta);
%     l_2 = @(theta) sin(2*pi*theta);
%     l_1_prime = @(theta) -2*pi*sin(theta*2*pi);
%     l_2_prime = @(theta) 2*pi*cos(theta*2*pi);
%     l_1_dprime = @(theta) -4*pi^2*cos(theta*2*pi);
%     l_2_dprime = @(theta) -4*pi^2*sin(theta*2*pi);
%     
    f_boundary = @(theta) l_1(theta).^2-l_2(theta).^2;

    K_boundary = @(theta_1, theta_2) -1/(2*pi) * ((l_1(theta_1)-l_1(theta_2)).*l_2_prime(theta_2)-(l_2(theta_1)-l_2(theta_2)).*l_1_prime(theta_2))./(sqrt(l_1_prime(theta_2).^2+l_2_prime(theta_2).^2).*((l_1(theta_1)-l_1(theta_2)).^2+(l_2(theta_1)-l_2(theta_2)).^2));
    K_boundary_same_point = @(theta) -1/(4*pi) * (l_2_prime(theta).*l_1_dprime(theta)-l_1_prime(theta).*l_2_dprime(theta))./(l_1_prime(theta).^2+l_2_prime(theta).^2).^(3/2);

    K_general = @(x, theta) -1/(2*pi) * ((x(1)-l_1(theta)).*l_2_prime(theta)-(x(2)-l_2(theta)).*l_1_prime(theta))./(sqrt(l_1_prime(theta).^2+l_2_prime(theta).^2).*((x(1)-l_1(theta)).^2+(x(2)-l_2(theta)).^2));

    p = 3;
    v = @(s) (1/p-1/2)*(1-2*s).^3+1/p*(2*s-1)+1/2;
    v_prime = @(s) 2/p-6*(1/p-1/2)*(1-2*s).^2;

    w = @(s) v(s).^p./(v(s).^p+v(1-s).^p);
    w_prime = @(s) p*((v_prime(s).*v(s).^(p-1))./(v(s).^p+v(1-s).^p)-(v(s).^(p-1).*v_prime(s)-v(1-s).^(p-1).*v_prime(1-s)).*v(s).^p./(v(s).^p+v(1-s).^p).^2);

    n = 80;
    ds = 1/n;
    s_mesh = linspace(0, 1-ds, n)';
    theta_mesh = w(s_mesh);

    figure;
    scatter(l_1(theta_mesh), l_2(theta_mesh));

    msk_diagonals = logical(diag(ones(n, 1)));

    [theta_2, theta_1] = meshgrid(theta_mesh);

    A = w_prime(s_mesh').*K_boundary(theta_1, theta_2).*sqrt(l_1_prime(theta_mesh').^2+l_2_prime(theta_mesh').^2)*ds; A(msk_diagonals) = 0;
    A(msk_diagonals) = w_prime(s_mesh).*K_boundary_same_point(theta_mesh).*sqrt(l_1_prime(theta_mesh).^2+l_2_prime(theta_mesh).^2)*ds + 1/2; 

    b = f_boundary(theta_mesh);
    phi_j = A\b;
    
%     figure;
% %     plot(A(1, 1:end))
%     hold on;
% %     plot(A(2, 1:end))
% %     plot(A(3, 1:end))
%     plot(A(10, 1:end))

    u_num = @(x, y)  transpose(w_prime(s_mesh).*K_general([x; y], theta_mesh).*sqrt(l_1_prime(theta_mesh).^2+l_2_prime(theta_mesh).^2).*phi_j) * ones(n, 1) * ds;
    
    % FFT coeff of phi computed
    fft_coeff_phi = fftshift(fft(phi_j))/n;
    
    R_fac = 2;
    n_R = n*R_fac;
    ds_R = ds/R_fac;
    s_R_mesh = linspace(0, 1-ds_R, n_R)';
    theta_R_mesh = w(s_R_mesh);
    
    padded_fft_coeffs = [zeros(floor((n_R-n)/2), 1); fft_coeff_phi; zeros(ceil((n_R-n)/2), 1)];
    phi_j_R = n_R*real(ifft(ifftshift(padded_fft_coeffs)));

    u_num_R = @(x, y)  transpose(w_prime(s__R_mesh).*K_general([x; y], theta_R_mesh).*sqrt(l_1_prime(theta_R_mesh).^2+l_2_prime(theta_R_mesh).^2).*phi_j_R) * ones(n_R, 1) * ds_R;
    
    if mod(n, 2) == 0
        freq_mesh = -n/2:(n/2-1);
    else
        freq_mesh = (-(n-1)/2):((n-1)/2);
    end
    
    msk_diagonals_R = logical(diag(ones(n_R, 1)));
    [theta_R_2, theta_R_1] = meshgrid(theta_R_mesh);
    A_R = w_prime(s_R_mesh').*K_boundary(theta_R_1, theta_R_2).*sqrt(l_1_prime(theta_R_mesh').^2+l_2_prime(theta_R_mesh').^2)*ds_R;
    
    A_R(msk_diagonals_R) = w_prime(s_R_mesh).*K_boundary_same_point(theta_R_mesh).*sqrt(l_1_prime(theta_R_mesh).^2+l_2_prime(theta_R_mesh).^2)*ds_R + 1/2;
    

    
    phi_j_R
    phi_j_exact = A_R\f_boundary(theta_R_mesh);
    phi_j_exact
    
    figure;
    plot(s_R_mesh, phi_j_exact)
    hold on;
    plot(s_R_mesh, phi_j_R)
    
    plot(s_mesh, phi_j)
%     legend('exact', 'estimate', 'og')
    err_Gamma = 1/n_R*sum(abs(phi_j-phi_j_exact(1:2:end)))
    
    figure;
    
        
    delta = 4*h;
    refined_theta_mesh = linspace(0, 1, n*100)';
    boundary_X = l_1(refined_theta_mesh);
    boundary_Y = l_2(refined_theta_mesh);

    R = R_cartesian_mesh_obj(x_start, x_end, y_start, y_end, h, boundary_X, boundary_Y);
    near_boundary_msk = compute_near_boundary_points(delta, R.R_X, R.R_Y, boundary_X, boundary_Y);

    interior_point_idxs = R.R_idxs(R.in_interior & ~near_boundary_msk);
    near_boundary_point_idxs = R.R_idxs(R.in_interior & near_boundary_msk);
    
    % numeric computation in interior
%     f_numeric = zeros(size(R.R_X));
%     for idx = interior_point_idxs'
%         f_numeric(idx) = transpose(w_prime(s_mesh).*K_general([R.R_X(idx); R.R_Y(idx)], theta_mesh).*sqrt(l_1_prime(theta_mesh).^2+l_2_prime(theta_mesh).^2).*phi_j) * ones(n, 1) * ds;
%     end
        
    
%     figure;
%     plot(s_R_mesh, phi_j_R, 'o')
%     hold on;
%     plot(s_mesh, phi_j, 'x')
%     plot(s_R_mesh, phi_j_R_exact)
%     legend('fourier', 'coarse', 'exact')

    err_(~R.in_interior) = nan;
end


function [near_boundary] = compute_near_boundary_points(delta, R_X, R_Y, boundary_X, boundary_Y) 
   %assumes uniform h
   h = R_X(1, 2) - R_X(1, 1);
   x_start = R_X(1, 1);
   y_start = R_Y(1, 1);
   n_x = size(R_X, 2);
   n_y = size(R_X, 1);
   
   near_boundary = false(size(R_X));
   
   straight_threshold = ceil(delta/h)+1;
   diagonal_threshold = ceil(delta/(sqrt(2)*h))+1;
   
   % finds closest cartesian point to each boundary point
   boundary_X_j = maintain_bounds(round((boundary_X-x_start)/h), 0, n_x-1);
   boundary_Y_j = maintain_bounds(round((boundary_Y-y_start)/h), 0, n_y-1);
   
   straight_left_right_idxs = maintain_bounds(boundary_X_j + (-straight_threshold:straight_threshold), 0, n_x-1);
   fixed_Y = boundary_Y_j + zeros(1, straight_threshold*2+1);
   near_boundary(sub2ind(size(near_boundary), fixed_Y(:)+1, straight_left_right_idxs(:)+1)) = true;
   
   straight_up_down_idxs = maintain_bounds(boundary_Y_j + (-straight_threshold:straight_threshold), 0, n_y-1);
   fixed_X = boundary_X_j + zeros(1, straight_threshold*2+1);
   near_boundary(sub2ind(size(near_boundary), straight_up_down_idxs(:)+1, fixed_X(:)+1)) = true;
   
   diag_1_X = maintain_bounds(boundary_X_j + (-diagonal_threshold:diagonal_threshold), 0, n_x-1);
   diag_1_Y = maintain_bounds(boundary_Y_j + (-diagonal_threshold:diagonal_threshold), 0, n_y-1);
   near_boundary(sub2ind(size(near_boundary), diag_1_Y(:)+1, diag_1_X(:)+1)) = true;
   
   diag_2_Y = maintain_bounds(boundary_Y_j + (diagonal_threshold:-1:-diagonal_threshold), 0, n_y-1);
   near_boundary(sub2ind(size(near_boundary), diag_2_Y(:)+1, diag_1_X(:)+1)) = true;
end

function idx_array = maintain_bounds(idx_array, min_bound, max_bound)
   idx_array(idx_array < min_bound) = min_bound;
   idx_array(idx_array > max_bound) = max_bound;
end