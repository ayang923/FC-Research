classdef C1_patch_obj < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        B
        B_c % B complement
        n_bound
        phi
        h
    end
    
    methods
        function obj = C1_patch_obj(M_p, J, n_bound, phi, f_B, f_B_c)
            obj.B = Q_patch_obj(M_p, J, n_bound, n_bound*2, 1/2, 1, 0, 1, f_B, phi);
            obj.B_c = Q_patch_obj(M_p, J, n_bound, n_bound, 0, 1/2, 1/2, 1, f_B_c, phi);
            
            obj.h = 1/(2*n_bound);
            
            obj.n_bound = n_bound;
            obj.phi = phi;
        end
        
        function [boundary_mesh_xi, boundary_mesh_eta] = boundary_mesh(obj)
            n_boundary = 10;
            boundary_mesh_xi = [zeros(n_boundary, 1); linspace(0, 1, n_boundary)'; ones(n_boundary, 1); linspace(1, 1/2, n_boundary)'; 1/2*ones(n_boundary, 1); linspace(1/2, 0, n_boundary)'];
            boundary_mesh_eta = [linspace(1/2, 1, n_boundary)'; ones(n_boundary, 1); linspace(1, 0, n_boundary)'; zeros(n_boundary, 1); linspace(0, 1/2, n_boundary)'; 1/2*ones(n_boundary, 1)];
        end
        
        function [boundary_mesh_x, boundary_mesh_y] = boundary_mesh_xy(obj)
            [boundary_mesh_xi, boundary_mesh_eta] = obj.boundary_mesh();
            [boundary_mesh_x, boundary_mesh_y] = obj.B.convert_to_XY(boundary_mesh_xi, boundary_mesh_eta);
        end
        
        function  [C1_fcont_patch_xi, C1_fcont_patch_eta_refined, C1_fcont_patch_eta_unrefined, B_fcont, B_c_minus_fcont] = FC(obj, C, n_r, d, A, Q, M, phi_B_normalization, phi_B_c_normalization)
            if ~isnan(phi_B_normalization)
                [XI, ETA] = obj.B.xi_eta_mesh;
                B_f_XY = obj.B.f_XY .* obj.B.phi(XI, ETA) ./ phi_B_normalization;
            else
                B_f_XY = obj.B.f_XY;
            end
            
            if ~isnan(phi_B_c_normalization)
                [XI, ETA] = obj.B_c.xi_eta_mesh;
                B_c_f_XY = obj.B_c.f_XY .* obj.B_c.phi(XI, ETA) ./ phi_B_c_normalization;
            else
                B_c_f_XY = obj.B_c.f_XY;
            end

            B_fcont = transpose(fcont_gram_blend_S(B_f_XY', d, A, Q));
            B_fc_unrefined = B_fcont(:, 1:n_r:n_r*C+1); % includes boundary point from where continuation is done
            
            C1_fcont_patch_xi = Q_patch_obj(obj.B.M_p, obj.B.J, C*n_r, obj.n_bound, obj.B.xi_start-C*obj.h, obj.B.xi_start, obj.B.eta_start, obj.B_c.eta_start, B_fcont(1:obj.n_bound+1, 1:C*n_r+1), obj.B.phi);
            
            if obj.n_bound < C
                B_c_minus_fcont = padarray(B_c_f_XY, [0 C - obj.n_bound], 0, 'pre') - B_fc_unrefined(obj.n_bound+1:end, :);
                B_c_minus_fcont_xi_start = obj.B_c.xi_start - (C-obj.n_bound)*obj.h;
            else
                [B_c_unrefined_f_XY, B_c_refined_f_XY] = obj.refine_B_c(B_c_f_XY, C, n_r, M);
                B_c_minus_fcont = B_c_refined_f_XY - B_fcont(obj.n_bound+1:end, 1:n_r*C+1);
                
                B_c_fcont_refined = fcont_gram_blend_S(B_c_minus_fcont, d, A, Q);
                B_c_fcont_unrefined = fcont_gram_blend_S(B_c_unrefined_f_XY, d, A, Q);
                
                C1_fcont_patch_eta_refined = Q_patch_obj(obj.B_c.M_p, obj.B_c.J, C*n_r, C*n_r, obj.B_c.xi_end-C*obj.h, obj.B_c.xi_end, obj.B_c.eta_start-C*obj.h, obj.B_c.eta_start, B_c_fcont_refined(1:C*n_r+1, :), obj.B_c.phi);
                C1_fcont_patch_eta_unrefined = Q_patch_obj(obj.B_c.M_p, obj.B_c.J, obj.n_bound-C, C*n_r, obj.B_c.xi_start, obj.B_c.xi_end-C*obj.h, obj.B_c.eta_start-C*obj.h, obj.B_c.eta_start, B_c_fcont_unrefined(1:C*n_r+1, :), obj.B_c.phi);
            end

%             B_c_fcont = fcont_gram_blend_S(B_c_minus_fcont, d, A, Q);
%             C1_fcont_patch_eta = Q_patch_obj(obj.B_c.M_p, obj.B_c.J, max(C, obj.n_bound), C*n_r, B_c_minus_fcont_xi_start, obj.B_c.xi_end, obj.B_c.eta_start-C*obj.h, obj.B_c.eta_start, B_c_fcont(1:C*n_r+1, :), obj.B_c.phi);
            
        end
        
        function [B_c_unrefined_f_XY, B_c_refined_f_XY] = refine_B_c(obj, B_c_f_XY, C, n_r, M)
            B_c_unrefined_f_XY = B_c_f_XY(:, 1:obj.n_bound-C+1);
            
            B_c_refined_xi_mesh = transpose((obj.B_c.xi_end - C*obj.h):(obj.h/n_r):obj.B_c.xi_end);
                                    
            B_c_refined_f_XY = zeros(obj.B_c.n_eta+1, length(B_c_refined_xi_mesh));
            half_M = floor((M+1)/2);
            for eta_j = 0:obj.n_bound
                for xi_j = obj.n_bound-C:(obj.n_bound-1)
                    if mod(M, 2) ~= 0
                        interpol_xi_j_mesh = transpose(xi_j-half_M+1:xi_j+half_M);
                    else
                        interpol_xi_j_mesh = transpose(xi_j-half_M:xi_j+half_M);
                    end

                    interpol_xi_j_mesh = shift_idx_mesh(interpol_xi_j_mesh, 0, obj.n_bound);
                    interpol_xi_mesh = obj.h*interpol_xi_j_mesh + obj.B_c.xi_start;
                    interpol_val = B_c_f_XY(eta_j+1, interpol_xi_j_mesh+1)';
                    
                    xi_eval_j = ((xi_j-(obj.n_bound-C))*n_r):((xi_j-(obj.n_bound-C)+1)*n_r);
                    
                    B_c_refined_f_XY(eta_j+1, xi_eval_j+1) = barylag([interpol_xi_mesh, interpol_val], B_c_refined_xi_mesh(xi_eval_j+1));
                end
            end
        end
        
        function [C1_B_norm, C1_B_c_norm, xi_norm, eta_norm] = compute_phi_normalization(obj, window_patch_xi, window_patch_eta)
            %computing corners of overlap on C2 patch
            C1_xi_corner = [max(convert_window_boundary(obj.B, window_patch_xi, window_patch_xi.eta_end, false, true)) ; max(convert_window_boundary(obj.B, window_patch_xi, window_patch_xi.xi_start, true, false))];
            C1_eta_corner = [max(convert_window_boundary(obj.B_c, window_patch_eta, window_patch_eta.xi_start, true, true)); min(convert_window_boundary(obj.B_c, window_patch_eta, window_patch_eta.eta_end, false, false))];
            
            % computation of phi
            C1_B_R_eta = C1_xi_corner(2) - obj.B.eta_start;
            C1_B_eta_0 = obj.B.eta_start;
            
            obj.B.phi = @(xi, eta) obj.B.phi_1D((eta-C1_B_eta_0)/C1_B_R_eta);
            
            C1_B_c_R_xi = C1_eta_corner(1) - obj.B_c.xi_start;
            C1_B_c_xi_0 = obj.B_c.xi_start;
            
            obj.B_c.phi = @(xi, eta) obj.B_c.phi_1D((xi - C1_B_c_xi_0)/C1_B_c_R_xi);
            
            % computing corners of overlap on window patches
            window_xi_corner = [max(convert_window_boundary(window_patch_xi, obj.B, obj.B.eta_start, false, true)); max(convert_window_boundary(window_patch_xi, obj.B, obj.B.xi_end, true, false))];
            window_eta_corner = [max(convert_window_boundary(window_patch_eta, obj.B_c, obj.B_c.xi_start, true, true)); max(convert_window_boundary(window_patch_eta, obj.B_c, obj.B_c.eta_end, false, false))];
            
            % computation of phi for window patches
            window_xi_R_xi = window_xi_corner(1) - window_patch_xi.xi_start;
            window_xi_xi_0 = window_patch_xi.xi_start;
            
            window_patch_xi.phi = @(xi, eta) window_patch_xi.phi_1D((xi-window_xi_xi_0)/window_xi_R_xi);

            window_eta_R_xi = window_eta_corner(1) - window_patch_eta.xi_start;
            window_eta_xi_0 = window_patch_eta.xi_start;

            window_patch_eta.phi = @(xi, eta) window_patch_eta.phi_1D((xi-window_eta_xi_0)/window_eta_R_xi) ;
            
            [C1_B_XI, C1_B_ETA] = obj.B.xi_eta_mesh();
            C1_B_norm = obj.B.phi(C1_B_XI, C1_B_ETA);
           
            [xi_overlap_XI, xi_overlap_ETA, xi_overlap_XI_j, xi_overlap_ETA_j] = compute_overlap_mesh(obj.B, C1_xi_corner, 3);
            [xi_overlap_X, xi_overlap_Y] = obj.B.convert_to_XY(xi_overlap_XI, xi_overlap_ETA);
            
            C1_B_norm = update_norm_arr(C1_B_norm, window_patch_xi, xi_overlap_X, xi_overlap_Y, xi_overlap_XI_j, xi_overlap_ETA_j, nan);

            [C1_B_c_XI, C1_B_c_ETA] = obj.B_c.xi_eta_mesh();
            C1_B_c_norm = obj.B_c.phi(C1_B_c_XI, C1_B_c_ETA);
            
            [eta_overlap_XI, eta_overlap_ETA, eta_overlap_XI_j, eta_overlap_ETA_j] = compute_overlap_mesh(obj.B_c, C1_eta_corner, 3);
            [eta_overlap_X, eta_overlap_Y] = obj.B_c.convert_to_XY(eta_overlap_XI, eta_overlap_ETA);
            
            C1_B_c_norm = update_norm_arr(C1_B_c_norm, window_patch_eta, eta_overlap_X, eta_overlap_Y, eta_overlap_XI_j, eta_overlap_ETA_j, nan);
            
            [xi_window_XI, xi_window_ETA] = window_patch_xi.xi_eta_mesh();
            xi_norm = window_patch_xi.phi(xi_window_XI, xi_window_ETA);
            
            [window_xi_overlap_XI, window_xi_overlap_ETA, window_xi_overlap_XI_j, window_xi_overlap_ETA_j] = compute_overlap_mesh(window_patch_xi, window_xi_corner, 3);
            [window_xi_overlap_X, window_xi_overlap_Y] = window_patch_xi.convert_to_XY(window_xi_overlap_XI, window_xi_overlap_ETA);
            
            xi_norm = update_norm_arr(xi_norm, obj.B, window_xi_overlap_X, window_xi_overlap_Y, window_xi_overlap_XI_j, window_xi_overlap_ETA_j, nan);
            
            [eta_window_XI, eta_window_ETA] = window_patch_eta.xi_eta_mesh();
            eta_norm = window_patch_eta.phi(eta_window_XI, eta_window_ETA);
                        
            [window_eta_overlap_XI, window_eta_overlap_ETA, window_eta_overlap_XI_j, window_eta_overlap_ETA_j] = compute_overlap_mesh(window_patch_eta, window_eta_corner, 3);
            [window_eta_overlap_X, window_eta_overlap_Y] = window_patch_eta.convert_to_XY(window_eta_overlap_XI, window_eta_overlap_ETA);
            
            eta_norm = update_norm_arr(eta_norm, obj.B_c, window_eta_overlap_X, window_eta_overlap_Y, window_eta_overlap_XI_j, window_eta_overlap_ETA_j, nan);
        end
    end
end

function idx_mesh = shift_idx_mesh(idx_mesh, min_bound, max_bound)
    if idx_mesh(1) < min_bound
            idx_mesh = idx_mesh + min_bound - idx_mesh(1);
    end
    if idx_mesh(end) > max_bound
            idx_mesh = idx_mesh + max_bound - idx_mesh(end);
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

