function [f_XY] = interpolate_cartesian(X, Y, X_prime, Y_prime, f_XY_prime, interior_msk)
%Interpolates (x_prime, y_prime, f_xy_prime) cartesian grid data on to the (x, y)
%non cartesian grid
% X and Y should be vectors of x-y values to interpolate on to and 
% X_prime, Y_prime should be 2D Cartesian matrices of x-y values.
% f_XY_prime is corresponding function values. interior_msk is the mask on
% the X_prime, Y_prime, and f_XY_prime for where data is known (data might
% not be perfect rectangle)
f_XY = zeros(size(X));
for i=1:length(X)
    x = X(i);
    y = Y(i);
    
    % computes closest j value to left
    x_j = floor((x+C*h_X)/h_X) + 1;
    y_j = floor((y+C*h_Y)/h_Y) + 1;
    
    % can check interior mask value at i, j to see if the interpolation is
    % within the bounds

    on_xi_grid = x-(-C*h_xi+(x_j-1)*h_xi) < eps;
    on_eta_grid = y-(-C*h_eta+(y_j-1)*h_eta) < eps;

    interpol_xi_j_mesh = transpose(x_j-(interpol_d+1)/2+1:x_j+(interpol_d+1)/2);
    interpol_eta_j_mesh = transpose(y_j-(interpol_d+1)/2+1:y_j+(interpol_d+1)/2);

    % boundary edge cases, shifts points to be within bounds
    interpol_xi_j_mesh = shift_idx_mesh(interpol_xi_j_mesh, 1, length(cont_xi_mesh));
    interpol_eta_j_mesh = shift_idx_mesh(interpol_eta_j_mesh, 1, length(cont_eta_mesh));

    interpol_xi_mesh = cont_xi_mesh(interpol_xi_j_mesh);
    interpol_eta_mesh = cont_eta_mesh(interpol_eta_j_mesh);
    % if corresponds to grid value, we don't need to interpolate
    if on_xi_grid && on_eta_grid
        patch_fcont(i) = fcont(x_j, y_j);
    else
        interpol_xi_points = zeros(size(interpol_xi_j_mesh));
        for vert_i=1:length(interpol_xi_points)
            interpol_coeff = polyfit(interpol_eta_mesh, fcont(interpol_xi_j_mesh(vert_i), interpol_eta_j_mesh), interpol_d);
            interpol_xi_points(vert_i) = polyval(interpol_coeff, y);
        end
        interpol_coeff = polyfit(interpol_xi_mesh, interpol_xi_points, interpol_d);
        patch_fcont(i) = polyval(interpol_coeff, x);
    end
end
end

