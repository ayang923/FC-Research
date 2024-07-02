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
    end
    
    methods
        function obj = Q_patch_obj(M_p, J, n_xi, n_eta, xi_start, xi_end, eta_start, eta_end, f_XY, phi)
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here
            obj.M_p = M_p;
            obj.J = J;
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
            h_xi =(obj.xi_end-obj.xi_start)./obj.n_xi;
            h_eta = (obj.eta_end-obj.eta_start) ./ obj.n_eta;
        end
        
        function mesh = xi_mesh(obj)
            % TODO: reduce roundoff error
            mesh = transpose(linspace(obj.xi_start, obj.xi_end, obj.n_xi+1));
        end
        
        function mesh = eta_mesh(obj)
            % TODO: reduce roundoff error
            mesh = transpose(linspace(obj.eta_start, obj.eta_end, obj.n_eta+1));
        end
        
        function [XI, ETA] = xi_eta_mesh(obj)
            [XI, ETA] = meshgrid(obj.xi_mesh(), obj.eta_mesh());
        end
        
        function [X, Y] = xy_mesh(obj)
            [XI, ETA] = obj.xi_eta_mesh();
            [X, Y] = obj.convert_to_XY(XI, ETA);
        end
        
        function [boundary_mesh_xi, boundary_mesh_eta] = boundary_mesh(obj)
            boundary_mesh_xi = [ones(obj.n_eta, 1)*obj.xi_start; obj.xi_mesh(); ones(obj.n_eta, 1)*obj.xi_end; flip(obj.xi_mesh())];
            boundary_mesh_eta = [obj.eta_mesh(); ones(obj.n_xi, 1)*obj.eta_end; flip(obj.eta_mesh); ones(obj.n_xi, 1)*obj.eta_start];
        end
        
        function [boundary_mesh_x, boundary_mesh_y] = boundary_mesh_xy(obj)
            [boundary_mesh_xi, boundary_mesh_eta] = obj.boundary_mesh();
            [boundary_mesh_x, boundary_mesh_y] = obj.convert_to_XY(boundary_mesh_xi, boundary_mesh_eta);
        end
        
        function [X, Y] = convert_to_XY(obj, XI, ETA)
            XY = obj.M_p(XI(:), ETA(:));
            X = reshape(XY(:, 1), size(XI));
            Y = reshape(XY(:, 2), size(ETA));
        end
        
        function patch_msk = in_patch(obj, xi, eta)
            eps = 0;
            patch_msk = xi >= obj.xi_start-eps & xi <= obj.xi_end+eps & eta >= obj.eta_start-eps & eta <= obj.eta_end+eps;
        end
        
        function patch_msk = in_patch_exact(obj, xi, eta)
            patch_msk =  xi >= obj.xi_start & xi <= obj.xi_end & eta >= obj.eta_start & eta <= obj.eta_end;
        end
        
        function patch_msk = in_patch_interior(obj, xi, eta)
            patch_msk =  xi > obj.xi_start & xi < obj.xi_end & eta > obj.eta_start & eta < obj.eta_end;
        end
        
        function [xi, eta] = round_boundary_points(obj, xi, eta)
            eps = 1e-14;
            
            xi(abs(xi - obj.xi_start) < eps) = obj.xi_start;
            xi(abs(xi-obj.xi_end) < eps) = obj.xi_end;
            eta(abs(eta - obj.eta_start) < eps) = obj.eta_start;
            eta(abs(eta-obj.eta_end) < eps) = obj.eta_end;
        end
        
        function [xi, eta, converged] = inverse_M_p(obj, x, y, initial_guesses)
            % need to change randomization, not giving  consistently good
            % results, probably ocnvergence issues.
            % initial_guess is [xi_1, xi_2, xi_3, ...; eta_1, eta_2, eta_3,
            % ...] matrix
            
            if isnan(initial_guesses)
                N = 20;
            
                N_segment = ceil(N/4);
                xi_mesh = transpose(linspace(obj.xi_start, obj.xi_end, N_segment+1));
                eta_mesh = transpose(linspace(obj.eta_start, obj.eta_end, N_segment+1));

                xi_initial = [xi_mesh(1:end-1); xi_mesh(1:end-1); ones(N_segment, 1)*obj.xi_start; ones(N_segment, 1)*obj.xi_end];
                eta_initial = [ones(N_segment, 1)*obj.eta_start; ones(N_segment, 1)*obj.eta_end; eta_mesh(1:end-1); eta_mesh(1:end-1)];
                
                initial_guesses = [xi_initial'; eta_initial'];                
            end

            err_guess = @(x, y, v) transpose(obj.M_p(v(1), v(2))) - [x; y];
            eps =  1e-14;
            if size(initial_guesses, 2) == 8
                disp('nan initial guess')
            end
            for k = 1:size(initial_guesses, 2)
                initial_guess = initial_guesses(:, k);
                [v_guess, converged] = newton_solve(@(v) err_guess(x, y, v), obj.J, initial_guess,eps, 100);
                
                [xi, eta] = obj.round_boundary_points(v_guess(1), v_guess(2));                                                
                if converged && obj.in_patch_exact(xi, eta)
                    return
                end
            end
        end
        
        function [f_xy, in_range] = locally_compute(obj, xi, eta, d, biperiodic)
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
            
            half_d =  floor((d+1)/2);
            if mod(d, 2) ~= 0
                interpol_xi_j_mesh = transpose(xi_j-half_d+1:xi_j+half_d);
                interpol_eta_j_mesh = transpose(eta_j-half_d+1:eta_j+half_d);
            else
                interpol_xi_j_mesh = transpose(xi_j-half_d:xi_j+half_d);
                interpol_eta_j_mesh = transpose(eta_j-half_d:eta_j+half_d);
            end

            if ~biperiodic
                interpol_xi_j_mesh = shift_idx_mesh(interpol_xi_j_mesh, 0, obj.n_xi);
                interpol_eta_j_mesh = shift_idx_mesh(interpol_eta_j_mesh, 0, obj.n_eta);
            end
            
            interpol_xi_mesh = h_xi*interpol_xi_j_mesh + obj.xi_start;
            interpol_eta_mesh = h_eta*interpol_eta_j_mesh + obj.eta_start;
            
            if biperiodic
                padded_f_XY = padarray(obj.f_XY, [d, d], 'both');
            end
            
            % first 1D interpolation
            interpol_eta_exact = zeros(d+1, 1);
            for vert_idx = 1:d+1
                if biperiodic
                    interpol_val = padded_f_XY(interpol_eta_j_mesh+d+1, interpol_xi_j_mesh(vert_idx)+d+1);
                else
                    interpol_val = obj.f_XY(interpol_eta_j_mesh+1, interpol_xi_j_mesh(vert_idx)+1);
                end
                interpol_eta_exact(vert_idx) = barylag([interpol_eta_mesh, interpol_val], eta);
            end
            
            % second 1D interpolation
            f_xy = barylag([interpol_xi_mesh, interpol_eta_exact], xi);
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

