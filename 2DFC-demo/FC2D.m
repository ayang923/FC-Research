function [R] = FC2D(f, h, curve_seq, eps_xi_eta, eps_xy, d, C, n_r, A, Q, M)
%FC2D Summary of this function goes here
%   Detailed explanation goes here
    patches = curve_seq.construct_patches(f, d, eps_xi_eta, eps_xy);
    
    x_min = curve_seq.first_curve.l_1(0);
    x_max = curve_seq.first_curve.l_1(0);
    y_min = curve_seq.first_curve.l_2(0);
    y_max = curve_seq.first_curve.l_2(0);
    
    FC_patches = cell(4*curve_seq.n_curves, 1);
    for i = 1:curve_seq.n_curves
        FC_patches{4*i-3} = patches{2*i-1}.FC(C, n_r, d, A, Q);
        [FC_patches{4*i-2}, FC_patches{4*i-1}, FC_patches{4*i}] = patches{2*i}.FC(C, n_r, d, A, Q, M);
        for j = 0:3
            x_min = min(FC_patches{4*i-j}.x_min, x_min);
            x_max = max(FC_patches{4*i-j}.x_max, x_max);
            y_min = min(FC_patches{4*i-j}.y_min, y_min);
            y_max = max(FC_patches{4*i-j}.y_max, y_max);
        end
    end

    [boundary_X, boundary_Y] = curve_seq.construct_boundary_mesh(n_r);
    R = R_cartesian_mesh_obj(x_min-h, x_max+h, y_min-h, y_max+h, h, boundary_X, boundary_Y);

    for i=1:length(FC_patches)
        R.interpolate_patch(FC_patches{i}, M);
    end

    R.fill_interior(f);
    
    R.compute_fc_coeffs();
    [R_X_err, R_Y_err, f_interpolation, interior_idx] = R.ifft_interpolation(R.h/2);
    f_exact = f(R_X_err, R_Y_err);
    disp(['fc error: ', num2str(max(abs(f_exact(interior_idx) - f_interpolation(interior_idx)), [], 'all'))])
end

