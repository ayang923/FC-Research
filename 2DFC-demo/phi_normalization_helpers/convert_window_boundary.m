function main_patch_boundary_mesh = convert_window_boundary(main_patch, window_patch, fixed_boundary_val, xi_window, xi_main)
    %Parameters: main_patch (convert coordinates to), window_patch (convert
    %coordinates from), fixed_boundary_val(0 or 1), which boundary we're
    %using in window_patch (xi, eta) coordinates, xi_window (true or false) if true, we are fixing xi in the window patch, if false, we are
    %fixing eta in the window patch. xi_main tells us which mesh we are
    %getting in the main coordinates
    if xi_window
        window_patch_eta_mesh = window_patch.eta_mesh();
        l = length(window_patch_eta_mesh);
        window_patch_xy_bounds = window_patch.M_p(fixed_boundary_val * ones(l, 1), window_patch_eta_mesh);
    else
        window_patch_xi_mesh = window_patch.xi_mesh();
        l = length(window_patch_xi_mesh);
        window_patch_xy_bounds = window_patch.M_p(window_patch_xi_mesh, fixed_boundary_val * ones(l, 1));
    end
    main_patch_boundary_coord = zeros(l, 2);
    initial_guesses = nan;
    for i = 1:l
        [main_patch_boundary_coord(i, 1), main_patch_boundary_coord(i, 2), converged] = main_patch.inverse_M_p(window_patch_xy_bounds(i, 1), window_patch_xy_bounds(i, 2), initial_guesses);        
        if ~converged
            warning("Nonconvergence in computing boundary mesh values")
        else
            initial_guesses = transpose(main_patch_boundary_coord(i, :));
        end

    end
    if xi_main
        main_patch_boundary_mesh = main_patch_boundary_coord(:, 1);
    else
        main_patch_boundary_mesh = main_patch_boundary_coord(:, 2);
    end
end