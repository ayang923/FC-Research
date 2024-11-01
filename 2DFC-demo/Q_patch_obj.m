classdef Q_patch_obj < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        M_p
        J
        n_xi
        n_eta
        xi_start
        xi_end
        eta_start
        eta_end
        f_XY
        x_min
        x_max
        y_min
        y_max
        phi_1D
        phi
        
        eps_xi_eta
        eps_xy
    end
    
    methods
        function obj = Q_patch_obj(M_p, J, eps_xi_eta, eps_xy, n_xi, n_eta, xi_start, xi_end, eta_start, eta_end, f_XY, phi)
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here
            obj.M_p = M_p;
            obj.J = J;
            
            obj.eps_xi_eta = eps_xi_eta;
            obj.eps_xy = eps_xy;
            
            obj.n_eta = n_eta;
            obj.n_xi = n_xi;
            obj.xi_start = xi_start;
            obj.xi_end = xi_end;
            obj.eta_start = eta_start;
            obj.eta_end = eta_end;
            obj.f_XY = f_XY;
            
            % sets XY bounds
            [XI, ETA] = obj.xi_eta_mesh();
            XY = M_p(XI(:), ETA(:));
            obj.x_min = min(XY(:, 1));
            obj.x_max = max(XY(:, 1));
            obj.y_min = min(XY(:, 2));
            obj.y_max = max(XY(:, 2));
            
            obj.phi_1D = @(x) erfc(6*(-2*x+1))/2;
            obj.phi = phi; %phi not defined until overlaps are given
        end

        function [h_xi, h_eta] = h_mesh(obj)
            h_xi =(obj.xi_end-obj.xi_start)./ (obj.n_xi-1);
            h_eta = (obj.eta_end-obj.eta_start) ./ (obj.n_eta-1);
        end
        
        function mesh = xi_mesh(obj)
            mesh = transpose(linspace(obj.xi_start, obj.xi_end, obj.n_xi));
        end
        
        function mesh = eta_mesh(obj)
            mesh = transpose(linspace(obj.eta_start, obj.eta_end, obj.n_eta));
        end
        
        function [XI, ETA] = xi_eta_mesh(obj)
            [XI, ETA] = meshgrid(obj.xi_mesh(), obj.eta_mesh());
        end
        
        function [X, Y] = xy_mesh(obj)
            [XI, ETA] = obj.xi_eta_mesh();
            [X, Y] = obj.convert_to_XY(XI, ETA);
        end
        
        function [boundary_mesh_xi, boundary_mesh_eta] = boundary_mesh(obj, pad_boundary)
            if pad_boundary
                [h_xi, h_eta] = obj.h_mesh;
            else
                h_xi = 0;
                h_eta = 0;
            end

            boundary_mesh_xi = [ones(obj.n_eta+1, 1)*(obj.xi_start-h_xi); obj.xi_mesh(); ones(obj.n_eta+1, 1)*(obj.xi_end+h_xi); flip(obj.xi_mesh()); obj.xi_start-h_xi];
            boundary_mesh_eta = [obj.eta_mesh(); ones(obj.n_xi+1, 1)*(obj.eta_end+h_eta); flip(obj.eta_mesh); ones(obj.n_xi+1, 1)*(obj.eta_start-h_eta); obj.eta_start];
        end
        
        function [boundary_mesh_x, boundary_mesh_y] = boundary_mesh_xy(obj, pad_boundary)
            [boundary_mesh_xi, boundary_mesh_eta] = obj.boundary_mesh(pad_boundary);
            [boundary_mesh_x, boundary_mesh_y] = obj.convert_to_XY(boundary_mesh_xi, boundary_mesh_eta);
        end
        
        function [X, Y] = convert_to_XY(obj, XI, ETA)
            XY = obj.M_p(XI(:), ETA(:));
            X = reshape(XY(:, 1), size(XI));
            Y = reshape(XY(:, 2), size(ETA));
        end
        
        function patch_msk = in_patch(obj, xi, eta)
            patch_msk = xi >= obj.xi_start & xi <= obj.xi_end & eta >= obj.eta_start & eta <= obj.eta_end;
        end
        
        function [xi, eta] = round_boundary_points(obj, xi, eta)
            xi(abs(xi - obj.xi_start) < obj.eps_xi_eta) = obj.xi_start;
            xi(abs(xi-obj.xi_end) < obj.eps_xi_eta) = obj.xi_end;
            eta(abs(eta - obj.eta_start) < obj.eps_xi_eta) = obj.eta_start;
            eta(abs(eta-obj.eta_end) < obj.eps_xi_eta) = obj.eta_end;
        end
        
        function [xi, eta, converged] = inverse_M_p(obj, x, y, initial_guesses)
            % initial_guess is [xi_1, xi_2, xi_3, ...; eta_1, eta_2, eta_3,
            % ...] matrix
            
            if isnan(initial_guesses)
                N = 20;
            
                N_segment = ceil(N/4);
                
                % TODO: this is just a boundary mesh
                xi_mesh = transpose(linspace(obj.xi_start, obj.xi_end, N_segment+1));
                eta_mesh = transpose(linspace(obj.eta_start, obj.eta_end, N_segment+1));

                xi_initial = [xi_mesh(1:end-1); xi_mesh(1:end-1); ones(N_segment, 1)*obj.xi_start; ones(N_segment, 1)*obj.xi_end];
                eta_initial = [ones(N_segment, 1)*obj.eta_start; ones(N_segment, 1)*obj.eta_end; eta_mesh(1:end-1); eta_mesh(1:end-1)];
                
                initial_guesses = [xi_initial'; eta_initial'];                
            end

            err_guess = @(x, y, v) transpose(obj.M_p(v(1), v(2))) - [x; y];

            for k = 1:size(initial_guesses, 2)
                initial_guess = initial_guesses(:, k);
                [v_guess, converged] = newton_solve(@(v) err_guess(x, y, v), obj.J, initial_guess, obj.eps_xy, 100);
                
                [xi, eta] = obj.round_boundary_points(v_guess(1), v_guess(2));                                                
                if converged && obj.in_patch(xi, eta)
                    return
                end
            end
        end

        function [f_xy, in_range] = locally_compute(obj, xi, eta, M)
            % check if xi, eta are within bounds of Q
            if ~obj.in_patch(xi, eta)
                f_xy = nan;
                in_range = false;
                return
            end
            
            in_range = true;
            [h_xi, h_eta] = obj.h_mesh();

            % j to the immediate left of point
            xi_j = floor((xi-obj.xi_start)/h_xi);
            eta_j = floor((eta-obj.eta_start)/h_eta);
            
            half_M =  floor(M/2);
            if mod(M, 2) ~= 0
                interpol_xi_j_mesh = transpose(xi_j-half_M:xi_j+half_M);
                interpol_eta_j_mesh = transpose(eta_j-half_M:eta_j+half_M);
            else
                interpol_xi_j_mesh = transpose(xi_j-half_M+1:xi_j+half_M);
                interpol_eta_j_mesh = transpose(eta_j-half_M+1:eta_j+half_M);
            end

            interpol_xi_j_mesh = shift_idx_mesh(interpol_xi_j_mesh, 0, obj.n_xi-1);
            interpol_eta_j_mesh = shift_idx_mesh(interpol_eta_j_mesh, 0, obj.n_eta-1);
            
            interpol_xi_mesh = h_xi*interpol_xi_j_mesh + obj.xi_start;
            interpol_eta_mesh = h_eta*interpol_eta_j_mesh + obj.eta_start;
            
            % first 1D interpolation
            interpol_xi_exact = zeros(M, 1);
            for horz_idx = 1:M
                mu = [mean(interpol_xi_mesh), std(interpol_xi_mesh)];
                interpol_val = obj.f_XY(interpol_eta_j_mesh(horz_idx)+1, interpol_xi_j_mesh+1)';
                interpol_xi_exact(horz_idx) = barylag([(interpol_xi_mesh-mu(1))/mu(2), interpol_val], (xi-mu(1))/mu(2));
            end
             % second 1D interpolation
            f_xy = barylag([interpol_eta_mesh, interpol_xi_exact], eta);
        end
        
        function [main_patch_phi, window_patch_phi] = compute_phi_normalization_xi_right(obj, window_patch)
            main_xi_corner = compute_xi_corner(obj, window_patch, true, window_patch.xi_start, true);
            window_xi_corner = compute_xi_corner(window_patch, obj, true, obj.xi_end, false);

            main_R_xi = main_xi_corner - obj.xi_end;
            main_xi_0 = obj.xi_end;
            obj.phi = @(xi, eta) obj.phi_1D((xi-main_xi_0)/main_R_xi);
            
            window_patch_phi = compute_phi_normalization_window(obj, window_patch, window_xi_corner);

            [XI_overlap, ETA_overlap, XI_j, ETA_j] = compute_xi_overlap_mesh(obj, main_xi_corner, true);
            [X_overlap, Y_overlap] = obj.convert_to_XY(XI_overlap, ETA_overlap);    
            
            [main_XI, main_ETA] = obj.xi_eta_mesh;
            main_patch_phi = obj.phi(main_XI, main_ETA);
            main_patch_phi = update_norm_arr(main_patch_phi, window_patch, X_overlap, Y_overlap, XI_j, ETA_j, nan);
        end
        
        function [main_patch_phi, window_patch_phi] = compute_phi_normalization_xi_left(obj, window_patch)
            main_xi_corner = compute_xi_corner(obj, window_patch, true, window_patch.xi_start, false);
            window_xi_corner = compute_xi_corner(window_patch, obj, true, obj.xi_start, false);

            main_R_xi = main_xi_corner - obj.xi_start;
            main_xi_0 = obj.xi_start;
            obj.phi = @(xi, eta) obj.phi_1D((xi-main_xi_0)/main_R_xi);
            
            window_patch_phi = compute_phi_normalization_window(obj, window_patch, window_xi_corner);

            [XI_overlap, ETA_overlap, XI_j, ETA_j] = compute_xi_overlap_mesh(obj, main_xi_corner, false);
            [X_overlap, Y_overlap] = obj.convert_to_XY(XI_overlap, ETA_overlap);                
            
            [main_XI, main_ETA] = obj.xi_eta_mesh;
            main_patch_phi = obj.phi(main_XI, main_ETA);
            main_patch_phi = update_norm_arr(main_patch_phi, window_patch, X_overlap, Y_overlap, XI_j, ETA_j, nan);        
        end
    
        function [main_patch_phi, window_patch_phi] = compute_phi_normalization_eta_up(obj, window_patch)
            main_eta_corner = compute_eta_corner(obj, window_patch, true, window_patch.xi_start, true);
            window_xi_corner = compute_xi_corner(window_patch, obj, false, obj.eta_end, false);

            main_R_eta = main_eta_corner - obj.eta_end;
            main_eta_0 = obj.eta_end;
            obj.phi = @(xi, eta) obj.phi_1D((eta-main_eta_0)/main_R_eta);
            
            window_patch_phi = compute_phi_normalization_window(obj, window_patch, window_xi_corner);

            [XI_overlap, ETA_overlap, XI_j, ETA_j] = compute_eta_overlap_mesh(obj, main_eta_corner, true);
            [X_overlap, Y_overlap] = obj.convert_to_XY(XI_overlap, ETA_overlap);         
            
            [main_XI, main_ETA] = obj.xi_eta_mesh;
            main_patch_phi = obj.phi(main_XI, main_ETA);
            main_patch_phi = update_norm_arr(main_patch_phi, window_patch, X_overlap, Y_overlap, XI_j, ETA_j, nan);
        end
        
        function [main_patch_phi, window_patch_phi] = compute_phi_normalization_eta_down(obj, window_patch)
            main_eta_corner = compute_eta_corner(obj, window_patch, true, window_patch.xi_start, false);
            window_xi_corner = compute_xi_corner(window_patch, obj, false, obj.eta_start, false);

            main_R_eta = main_eta_corner - obj.eta_start;
            main_eta_0 = obj.eta_start;
            obj.phi = @(xi, eta) obj.phi_1D((eta-main_eta_0)/main_R_eta);
            
            window_patch_phi = compute_phi_normalization_window(obj, window_patch, window_xi_corner);

            [XI_overlap, ETA_overlap, XI_j, ETA_j] = compute_eta_overlap_mesh(obj, main_eta_corner, true);
            [X_overlap, Y_overlap] = obj.convert_to_XY(XI_overlap, ETA_overlap);        
            
            [main_XI, main_ETA] = obj.xi_eta_mesh;
            main_patch_phi = obj.phi(main_XI, main_ETA);
            main_patch_phi = update_norm_arr(main_patch_phi, window_patch, X_overlap, Y_overlap, XI_j, ETA_j, nan);
        end
    end
end

function [main_xi_corner] = compute_xi_corner(main_patch, window_patch, window_fix_xi, window_fixed_edge, window_patch_right)
%COMPUTE_XI_CORNER Summary of this function goes here
%   Detailed explanation goes here
        [h_xi_window, h_eta_window] = window_patch.h_mesh;
        
        if window_fix_xi
            window_xi_edge = window_fixed_edge;
            window_eta_edge = window_patch.eta_start;
        else
            window_xi_edge = window_patch.xi_start;
            window_eta_edge = window_fixed_edge;
        end
                
        if window_patch_right
            main_xi_corner = main_patch.xi_end;
        else
            main_xi_corner = main_patch.xi_start;
        end

        while true
            window_xy_edge = window_patch.M_p(window_xi_edge, window_eta_edge);
            [main_xi_edge, main_eta_edge, converged] = main_patch.inverse_M_p(window_xy_edge(1, 1), window_xy_edge(1, 2), nan);
            if ~converged
                 warning("Nonconvergence in computing boundary mesh values");
                 break;
            end
            if main_eta_edge > main_patch.eta_end || main_eta_edge < main_patch.eta_start
                break;
            end
            if (main_xi_edge < main_xi_corner && window_patch_right) || (main_xi_edge > main_xi_corner && ~window_patch_right)
                main_xi_corner = main_xi_edge;
            end
            if window_fix_xi
                window_eta_edge = window_eta_edge + h_eta_window;
            else
                window_xi_edge = window_xi_edge + h_xi_window;
            end
        end
end

function [XI_overlap, ETA_overlap, XI_j, ETA_j] = compute_xi_overlap_mesh(main_patch, xi_corner, window_patch_right)
%COMPUTE_XI_OVERLAP_MESH Summary of this function goes here
%   Detailed explanation goes here
    [h_xi, h_eta] = main_patch.h_mesh;
    
    if window_patch_right
        xi_corner_j = ceil((xi_corner - main_patch.xi_start) / h_xi);
        [XI_j, ETA_j] = meshgrid(xi_corner_j:(main_patch.n_xi-1), 0:(main_patch.n_eta-1));        
    else
        xi_corner_j = floor((xi_corner - main_patch.xi_start) / h_xi);
        [XI_j, ETA_j] = meshgrid(0:xi_corner_j, 0:(main_patch.n_eta-1));
    end
    
    XI_overlap = XI_j * h_xi + main_patch.xi_start;
    ETA_overlap = ETA_j * h_eta + main_patch.eta_start;
end

function [main_eta_corner] = compute_eta_corner(main_patch, window_patch, window_fix_xi, window_fixed_edge, window_patch_up)
%COMPUTE_XI_CORNER Summary of this function goes here
%   Detailed explanation goes here
        [h_xi_window, h_eta_window] = window_patch.h_mesh;
        
        if window_fix_xi
            window_xi_edge = window_fixed_edge;
            window_eta_edge = window_patch.eta_start;
        else
            window_xi_edge = window_patch.xi_start;
            window_eta_edge = window_fixed_edge;
        end
                
        if window_patch_up
            main_eta_corner = main_patch.eta_end;
        else
            main_eta_corner = main_patch.eta_start;
        end

        while true
            window_xy_edge = window_patch.M_p(window_xi_edge, window_eta_edge);
            [main_xi_edge, main_eta_edge, converged] = main_patch.inverse_M_p(window_xy_edge(1, 1), window_xy_edge(1, 2), nan);
            if ~converged
                 warning("Nonconvergence in computing boundary mesh values");
                 break;
            end
            if main_xi_edge > main_patch.xi_end || main_xi_edge < main_patch.xi_start
                break;
            end
            if (main_eta_edge < main_eta_corner && window_patch_up) || (main_eta_edge > main_eta_corner && ~window_patch_up)
                main_eta_corner = main_eta_edge;
            end
            if window_fix_xi
                window_eta_edge = window_eta_edge + h_eta_window;
            else
                window_xi_edge = window_xi_edge + h_xi_window;
            end
        end
end



function [XI_overlap, ETA_overlap, XI_j, ETA_j] = compute_eta_overlap_mesh(main_patch, eta_corner, window_patch_up)
%COMPUTE_XI_OVERLAP_MESH Summary of this function goes here
%   Detailed explanation goes here
    [h_xi, h_eta] = main_patch.h_mesh;
    
    if window_patch_up
        eta_corner_j = ceil((eta_corner - main_patch.eta_start) / h_eta);
        [XI_j, ETA_j] = meshgrid(0:(main_patch.n_xi-1), eta_corner_j:(main_patch.n_eta-1));        
    else
        eta_corner_j = floor((eta_corner - main_patch.eta_start) / h_eta);
        [XI_j, ETA_j] = meshgrid(0:(main_patch.n_xi-1), 0:eta_corner_j);
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
            
            if converged
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

function window_patch_phi = compute_phi_normalization_window(main_patch, window_patch, window_xi_corner)
    window_R_xi = window_xi_corner - window_patch.xi_start;
    window_xi_0 = window_patch.xi_start;
    window_patch.phi = @(xi, eta) window_patch.phi_1D((xi-window_xi_0)/window_R_xi);

    [XI_overlap_window, ETA_overlap_window, XI_j_window, ETA_j_window] = compute_xi_overlap_mesh(window_patch, window_xi_corner, false);
    [X_overlap_window, Y_overlap_window] = window_patch.convert_to_XY(XI_overlap_window, ETA_overlap_window);
    
    [window_XI, window_ETA] = window_patch.xi_eta_mesh;
    window_patch_phi = window_patch.phi(window_XI, window_ETA);
    window_patch_phi = update_norm_arr(window_patch_phi, main_patch, X_overlap_window, Y_overlap_window, XI_j_window, ETA_j_window, nan);
end




