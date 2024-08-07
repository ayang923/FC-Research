classdef C2_patch_old_obj < Q_patch_obj
    %C2_PATCH_OBJ Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
        function obj = C2_patch_old_obj(M_p, J, n_xi, n_eta, xi_start, xi_end, eta_start, eta_end, f_XY, phi)
            %C2_PATCH_OBJ Construct an instance of this class
            %   Detailed explanation goes here
            obj = obj@Q_patch_obj(M_p, J, n_xi, n_eta, xi_start, xi_end, eta_start, eta_end, f_XY, phi);
        end
                
        function C2_fcont_patch = FC(obj, C, d, A, Q, phi_normalization)
            h_xi = (obj.xi_end-obj.xi_start)/obj.n_xi;
            h_eta = (obj.eta_end-obj.eta_start)/obj.n_eta;
            
            [XI, ETA] = obj.xi_eta_mesh();
            
            if ~isnan(phi_normalization)
                fcont = fcont_gram_blend_C2(obj.f_XY.*obj.phi(XI, ETA)./phi_normalization, d, A, Q);
            else
                fcont = fcont_gram_blend_C2(obj.f_XY, d, A, Q);
            end
            
            C2_fcont_patch = C2_patch_obj(obj.M_p, obj.J, C+obj.n_xi, C+obj.n_eta, obj.xi_start-(C)*h_xi, obj.xi_end, obj.eta_start-(C)*h_eta, obj.eta_end, fcont, obj.phi);
        end
        
        function [C2_norm, xi_norm, eta_norm] = compute_phi_normalization(obj, window_patch_xi, window_patch_eta)
            %computing corners of overlap on C2 patch
            C2_xi_corner = [min(convert_window_boundary(obj, window_patch_xi, window_patch_xi.xi_start, true, true)) ; max(convert_window_boundary(obj, window_patch_xi, window_patch_xi.eta_end, false, false))];
            C2_eta_corner = [max(convert_window_boundary(obj, window_patch_eta, window_patch_eta.eta_end, false, true)); min(convert_window_boundary(obj, window_patch_eta, window_patch_eta.xi_start, true, false))];
            
            % computation of phi
            C2_R_xi = C2_xi_corner(1) - obj.xi_end;
            C2_R_eta = C2_eta_corner(2) - obj.eta_end;
            C2_xi_0 = obj.xi_end;
            C2_eta_0 = obj.eta_end;
            
            obj.phi = @(xi, eta) obj.phi_1D((xi-C2_xi_0)/C2_R_xi) .* obj.phi_1D((eta-C2_eta_0)./C2_R_eta);            
            
            % computing corners of overlap on window patches
            window_xi_corner = [max(convert_window_boundary(window_patch_xi, obj, obj.xi_end, true, true)); max(convert_window_boundary(window_patch_xi, obj, obj.eta_end, false, false))];
            window_eta_corner = [max(convert_window_boundary(window_patch_eta, obj, obj.eta_end, false, true)); max(convert_window_boundary(window_patch_eta, obj, obj.xi_end, true, false))];
            
            % computation of phi for window patches
            window_xi_R_xi = window_xi_corner(1) - window_patch_xi.xi_start;
            window_xi_xi_0 = window_patch_xi.xi_start;
            
            window_patch_xi.phi = @(xi, eta) obj.phi_1D((xi-window_xi_xi_0)/window_xi_R_xi);
            
            window_eta_R_xi = window_eta_corner(1) - window_patch_eta.xi_start;
            window_eta_xi_0 = window_patch_eta.xi_start;
            
            window_patch_eta.phi = @(xi, eta) obj.phi_1D((xi-window_eta_xi_0)/window_eta_R_xi) ;
            
            [C2_XI, C2_ETA] = obj.xi_eta_mesh();
            C2_norm = obj.phi(C2_XI, C2_ETA);
           
            [xi_overlap_XI, xi_overlap_ETA, xi_overlap_XI_j, xi_overlap_ETA_j] = compute_overlap_mesh(obj, C2_xi_corner, 4);
            [xi_overlap_X, xi_overlap_Y] = obj.convert_to_XY(xi_overlap_XI, xi_overlap_ETA);            
            
            [eta_overlap_XI, eta_overlap_ETA, eta_overlap_XI_j, eta_overlap_ETA_j] = compute_overlap_mesh(obj, C2_eta_corner, 2);
            [eta_overlap_X, eta_overlap_Y] = obj.convert_to_XY(eta_overlap_XI, eta_overlap_ETA);          
            
            C2_norm = update_norm_arr(C2_norm, window_patch_xi, xi_overlap_X, xi_overlap_Y, xi_overlap_XI_j, xi_overlap_ETA_j, nan);
            C2_norm = update_norm_arr(C2_norm, window_patch_eta, eta_overlap_X, eta_overlap_Y, eta_overlap_XI_j, eta_overlap_ETA_j, nan);

            [xi_window_XI, xi_window_ETA] = window_patch_xi.xi_eta_mesh();
            xi_norm = window_patch_xi.phi(xi_window_XI, xi_window_ETA);
            
            [window_xi_overlap_XI, window_xi_overlap_ETA, window_xi_overlap_XI_j, window_xi_overlap_ETA_j] = compute_overlap_mesh(window_patch_xi, window_xi_corner, 3);
            [window_xi_overlap_X, window_xi_overlap_Y] = window_patch_xi.convert_to_XY(window_xi_overlap_XI, window_xi_overlap_ETA);
            
            figure;
            [X, Y] = obj.xy_mesh;
            scatter(X(:), Y(:));
            hold on;
            [X, Y] = window_patch_xi.xy_mesh;
            scatter(X(:), Y(:));
            scatter(window_xi_overlap_X(:), window_xi_overlap_Y(:));
            legend('C2', 'window', 'overlap')
            
            xi_norm = update_norm_arr(xi_norm, obj, window_xi_overlap_X, window_xi_overlap_Y, window_xi_overlap_XI_j, window_xi_overlap_ETA_j, nan);
            
            [eta_window_XI, eta_window_ETA] = window_patch_eta.xi_eta_mesh();
            eta_norm = window_patch_eta.phi(eta_window_XI, eta_window_ETA);
                        
            [window_eta_overlap_XI, window_eta_overlap_ETA, window_eta_overlap_XI_j, window_eta_overlap_ETA_j] = compute_overlap_mesh(window_patch_eta, window_eta_corner, 3);
            [window_eta_overlap_X, window_eta_overlap_Y] = window_patch_eta.convert_to_XY(window_eta_overlap_XI, window_eta_overlap_ETA);
            
            eta_norm = update_norm_arr(eta_norm, obj, window_eta_overlap_X, window_eta_overlap_Y, window_eta_overlap_XI_j, window_eta_overlap_ETA_j, nan);
        end
    end
end

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

function [XI_overlap, ETA_overlap, XI_j, ETA_j] = compute_overlap_mesh(main_patch, xi_eta_corner, quadrant_overlap)
    [h_xi, h_eta] = main_patch.h_mesh();
    j_corner =  (xi_eta_corner - [main_patch.xi_start; main_patch.eta_start])./ [h_xi; h_eta];
    
    if quadrant_overlap == 1
        j_corner_grid = [ceil(j_corner(1)); ceil(j_corner(2))];
        [XI_j, ETA_j] = meshgrid(max([j_corner_grid(1); 0]):main_patch.n_xi, max([j_corner_grid(2); 0]):main_patch.n_eta);
    elseif quadrant_overlap == 2
        j_corner_grid = [floor(j_corner(1)); ceil(j_corner(2))];
        [XI_j, ETA_j] = meshgrid(0:min([main_patch.n_xi; j_corner_grid(1)]), max([j_corner_grid(2); 0]):main_patch.n_eta);
    elseif quadrant_overlap == 3
        j_corner_grid = [floor(j_corner(1)); floor(j_corner(2))];
        [XI_j, ETA_j] = meshgrid(0:min([main_patch.n_xi; j_corner_grid(1)]), 0:min([main_patch.n_eta; j_corner_grid(2)]));
    elseif quadrant_overlap == 4
        j_corner_grid = [ceil(j_corner(1)); floor(j_corner(2))];
        [XI_j, ETA_j] = meshgrid(max([j_corner_grid(1); 0]):main_patch.n_xi, 0:min([main_patch.n_eta; j_corner_grid(2)]));
    else
        error("Invalid Quadrant Number")
    end
    XI_overlap = XI_j * h_xi + main_patch.xi_start;
    ETA_overlap = ETA_j * h_eta + main_patch.eta_start;
end

function [norm_arr] = update_norm_arr(norm_arr, window_patch, overlap_X, overlap_Y, overlap_XI_j, overlap_ETA_j, initial_guesses)
    for i = 1:size(overlap_X, 1)
        if mod(i, 2) == 1
            j_lst = 1:size(overlap_X, 2);
        else
            j_lst = size(overlap_X, 2):-1:1;
        end
        
        for j = j_lst
            [window_patch_xi, window_patch_eta, converged] = window_patch.inverse_M_p(overlap_X(i, j), overlap_Y(i, j), initial_guesses);
            
            if converged && window_patch.in_patch(window_patch_xi, window_patch_eta)
                xi_i = overlap_XI_j(i, j) + 1;
                xi_j = overlap_ETA_j(i, j) + 1;
                norm_arr(xi_j, xi_i) = norm_arr(xi_j, xi_i) + window_patch.phi(window_patch_xi, window_patch_eta);
            elseif ~converged
                warning("Nonconvergence in computing C2_norm")
            end
            
            if converged
                initial_guesses = [window_patch_xi; window_patch_eta]; % using previous point as next initial guess
            end
        end
    end
end