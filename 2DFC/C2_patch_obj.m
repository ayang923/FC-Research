classdef C2_patch_obj < Q_patch_obj
    %C2_PATCH_OBJ Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
        function obj = C2_patch_obj(M_p, J, n_xi, n_eta, xi_start, xi_end, eta_start, eta_end, f_XY)
            %C2_PATCH_OBJ Construct an instance of this class
            %   Detailed explanation goes here
            obj = obj@Q_patch_obj(M_p, J, n_xi, n_eta, xi_start, xi_end, eta_start, eta_end, f_XY);
        end
        
        function phi_fh = phi(obj, xi, eta)
            xi_0 = obj.xi_start;
            eta_0 = obj.eta_start;
            R_xi = (obj.xi_end-obj.xi_start);
            R_eta = (obj.eta_end-obj.eta_start);

            in_rectangle =obj.in_patch(xi, eta);
            phi_fh = zeros(size(xi));
            
            phi_fh(in_rectangle) = exp(-1./(1-(1/R_xi.^2).*(xi(in_rectangle)-xi_0).^2)).*exp(-1./(1-(1/R_eta.^2).*(eta(in_rectangle)-eta_0).^2));
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
            
            C2_fcont_patch = C2_patch_obj(obj.M_p, obj.J, C+obj.n_xi, C+obj.n_eta, obj.xi_start-C*h_xi, obj.xi_end, obj.eta_start-C*h_eta, obj.eta_end, fcont);
        end
        
        function [C2_norm, xi_norm, eta_norm] = compute_phi_normalization(obj, window_patch_xi, window_patch_eta)
            window_patch_xi_eta_mesh = window_patch_xi.eta_mesh();
            window_patch_eta_eta_mesh = window_patch_eta.eta_mesh();
            
            l_xi = 9; %length(window_patch_xi_eta_mesh);
            l_eta =9; % length(window_patch_eta_eta_mesh);
            
            window_patch_xi_xy_bound = window_patch_xi.M_p(zeros(l_xi, 1), window_patch_xi_eta_mesh(1:l_xi));
            window_patch_eta_xy_bound = window_patch_eta.M_p(zeros(l_eta, 1), window_patch_eta_eta_mesh(1:l_eta));
           
            C2_xi_bounds = zeros(l_xi, 2);
            C2_eta_bounds = zeros(l_eta, 2);
            for i=1:l_xi
                [C2_xi_bounds(i, 1), C2_xi_bounds(i, 2), converged] = obj.inverse_M_p(window_patch_xi_xy_bound(i, 1), window_patch_xi_xy_bound(i, 2));
                if ~converged
                    disp("not converged")
                end
            end
            for i=1:l_eta
                [C2_eta_bounds(i, 1), C2_eta_bounds(i, 2), converged] = obj.inverse_M_p(window_patch_eta_xy_bound(i, 1), window_patch_eta_xy_bound(i, 2));
                if ~converged
                    disp("not converged")
                end
            end
       
            C2_xi_max = [min(C2_xi_bounds(:, 1)) ; max(C2_xi_bounds(:, 2))];
            C2_eta_max = [max(C2_eta_bounds(:, 1)); min(C2_eta_bounds(:, 2))];
            
            [C2_XI, C2_ETA] = obj.xi_eta_mesh();
            C2_norm = obj.phi(C2_XI, C2_ETA);
            
            [C2_h_xi, C2_h_eta] = obj.h_mesh();
            
            C2_xi_j_bounds = (C2_xi_max - [obj.xi_start; obj.eta_start])./ [C2_h_xi; C2_h_eta];
            C2_xi_j_bounds = [max([ceil(C2_xi_j_bounds(1)); 0]); min([floor(C2_xi_j_bounds(2)); obj.n_eta])];
            
            C2_eta_j_bounds = (C2_eta_max - [obj.xi_start; obj.eta_start])./ [C2_h_xi; C2_h_eta];
            C2_eta_j_bounds = [min([floor(C2_eta_j_bounds(1)); obj.n_xi]); max([ceil(C2_eta_j_bounds(2)); 0])];
            
            [xi_overlap_XI_j, xi_overlap_ETA_j] = meshgrid(C2_xi_j_bounds(1):(obj.n_xi), 0:(C2_xi_j_bounds(2)));
            xi_overlap_XI = xi_overlap_XI_j * C2_h_xi + obj.xi_start;
            xi_overlap_ETA = xi_overlap_ETA_j * C2_h_eta + obj.eta_start;
            xi_overlap_XY = obj.M_p(xi_overlap_XI(:), xi_overlap_ETA(:));
            xi_overlap_X = reshape(xi_overlap_XY(:, 1), size(xi_overlap_XI));
            xi_overlap_Y = reshape(xi_overlap_XY(:, 2), size(xi_overlap_XI));
            
            [eta_overlap_XI_j, eta_overlap_ETA_j] = meshgrid(0:C2_eta_j_bounds(1), C2_eta_j_bounds(2):obj.n_eta);
            eta_overlap_XI = eta_overlap_XI_j * C2_h_xi + obj.xi_start;
            eta_overlap_ETA = eta_overlap_ETA_j * C2_h_eta + obj.eta_start;
            eta_overlap_XY = obj.M_p(eta_overlap_XI(:), eta_overlap_ETA(:));
            eta_overlap_X = reshape(eta_overlap_XY(:, 1), size(eta_overlap_XI));
            eta_overlap_Y = reshape(eta_overlap_XY(:, 2), size(eta_overlap_XI));
            
            max_window_xi = [window_patch_xi.xi_start; window_patch_xi.eta_start];
            for i = 1:size(xi_overlap_XI, 1)
                for j = 1:size(xi_overlap_XI, 2)
                    C2_xi_i = xi_overlap_XI_j(i, j) + 1;
                    C2_xi_j = xi_overlap_ETA_j(i, j) + 1;
                    [window_patch_xi_xi, window_patch_xi_eta, converged] = window_patch_xi.inverse_M_p(xi_overlap_X(i, j), xi_overlap_Y(i, j));
                    if ~converged
                        disp("not converged C2")
                    end
                    if converged && window_patch_xi.in_patch(window_patch_xi_xi, window_patch_xi_eta)
                        C2_norm(C2_xi_j, C2_xi_i) = C2_norm(C2_xi_j, C2_xi_i) + window_patch_xi.window_phi(window_patch_xi_xi, window_patch_xi_eta);

                        if window_patch_xi_xi > max_window_xi(1) && window_patch_xi_xi <= window_patch_xi.xi_end
                            max_window_xi(1) = window_patch_xi_xi;
                        end
                        if window_patch_xi_eta > max_window_xi(2) && window_patch_xi_eta <= window_patch_xi.eta_end
                            max_window_xi(2) = window_patch_xi_eta;
                        end
                    end
                end
            end
            
            max_window_eta = [window_patch_eta.xi_start; window_patch_eta.eta_start];
            for i = 1:size(eta_overlap_XI, 1)
                for j = 1:size(eta_overlap_XI, 2)
                    C2_eta_i = eta_overlap_XI_j(i, j) + 1;
                    C2_eta_j = eta_overlap_ETA_j(i, j) + 1;
                    [window_patch_eta_xi, window_patch_eta_eta, converged] = window_patch_eta.inverse_M_p(eta_overlap_X(i, j), eta_overlap_Y(i, j));
                    if ~converged
                        disp("not converged C2")
                    end
                    if converged && window_patch_eta.in_patch(window_patch_eta_xi, window_patch_eta_eta)
                        C2_norm(C2_eta_j, C2_eta_i) = C2_norm(C2_eta_j, C2_eta_i) + window_patch_eta.window_phi(window_patch_eta_xi, window_patch_eta_eta);

                        if window_patch_eta_xi > max_window_eta(1) && window_patch_eta_xi <= window_patch_eta.xi_end
                            max_window_eta(1) = window_patch_eta_xi;
                        end
                        if window_patch_eta_eta > max_window_eta(2) && window_patch_eta_eta <= window_patch_eta.eta_end
                            max_window_eta(2) = window_patch_eta_eta;
                        end
                    end
                end
            end
            
            [xi_window_XI, xi_window_ETA] = window_patch_xi.xi_eta_mesh();
            xi_norm = window_patch_xi.window_phi(xi_window_XI, xi_window_ETA);
            
            [window_xi_h_xi, window_xi_h_eta] = window_patch_xi.h_mesh();
            window_xi_j_bounds = (max_window_xi - [window_patch_xi.xi_start; window_patch_xi.eta_start])./ [window_xi_h_xi; window_xi_h_eta];
            window_xi_j_bounds = [min([ceil(window_xi_j_bounds(1)); obj.n_xi]); min([ceil(window_xi_j_bounds(2)); obj.n_eta])];
            
            [window_xi_overlap_XI_j, window_xi_overlap_ETA_j] = meshgrid(0:window_xi_j_bounds(1), 0:window_xi_j_bounds(2));
            window_xi_overlap_XI = window_xi_overlap_XI_j * window_xi_h_xi + window_patch_xi.xi_start;
            window_xi_overlap_ETA = window_xi_overlap_ETA_j * window_xi_h_eta + window_patch_xi.eta_start;
            window_xi_overlap_XY = window_patch_xi.M_p(window_xi_overlap_XI(:), window_xi_overlap_ETA(:));
            window_xi_overlap_X = reshape(window_xi_overlap_XY(:, 1), size(window_xi_overlap_XI));
            window_xi_overlap_Y = reshape(window_xi_overlap_XY(:, 2), size(window_xi_overlap_XI));
            
            for i = 1:size(window_xi_overlap_XI, 1)
                for j = 1:size(window_xi_overlap_XI, 2)
                    window_xi_i = window_xi_overlap_XI_j(i, j) + 1;
                    window_xi_j = window_xi_overlap_ETA_j(i, j) + 1;
                    [C2_xi, C2_eta, converged] = obj.inverse_M_p(window_xi_overlap_X(i, j), window_xi_overlap_Y(i, j));
                    if ~converged
                        disp("not converged xi")
                    end
                    if converged && obj.in_patch(C2_xi, C2_eta)
                        xi_norm(window_xi_j, window_xi_i) =  xi_norm(window_xi_j, window_xi_i) + obj.phi(C2_xi, C2_eta);
                    end
                end
            end
            
            [eta_window_XI, eta_window_ETA] = window_patch_eta.xi_eta_mesh();
            eta_norm = window_patch_eta.window_phi(eta_window_XI, eta_window_ETA);
            
            [window_eta_h_xi, window_eta_h_eta] = window_patch_eta.h_mesh();
            window_eta_j_bounds = (max_window_eta - [window_patch_eta.xi_start; window_patch_eta.eta_start])./ [window_eta_h_xi; window_eta_h_eta];
            window_eta_j_bounds = [min([ceil(window_eta_j_bounds(1)); obj.n_xi]); min([ceil(window_eta_j_bounds(2)); obj.n_eta])];
            
            [window_eta_overlap_XI_j, window_eta_overlap_ETA_j] = meshgrid(0:window_eta_j_bounds(1), 0:window_eta_j_bounds(2));
            window_eta_overlap_XI = window_eta_overlap_XI_j * window_eta_h_xi + window_patch_eta.xi_start;
            window_eta_overlap_ETA = window_eta_overlap_ETA_j * window_eta_h_eta + window_patch_eta.eta_start;
            window_eta_overlap_XY = window_patch_eta.M_p(window_eta_overlap_XI(:), window_eta_overlap_ETA(:));
            window_eta_overlap_X = reshape(window_eta_overlap_XY(:, 1), size(window_eta_overlap_XI));
            window_eta_overlap_Y = reshape(window_eta_overlap_XY(:, 2), size(window_eta_overlap_XI));
            
            for i = 1:size(window_eta_overlap_XI, 1)
                for j = 1:size(window_eta_overlap_XI, 2)
                    window_eta_i = window_eta_overlap_XI_j(i, j) + 1;
                    window_eta_j = window_eta_overlap_ETA_j(i, j) + 1;
                    [C2_xi, C2_eta, converged] = obj.inverse_M_p(window_eta_overlap_X(i, j), window_eta_overlap_Y(i, j));
                    if ~converged
                        disp("not converged eta")
                    end
                    if converged && obj.in_patch(C2_xi, C2_eta)
                        eta_norm(window_eta_j, window_eta_i) =  eta_norm(window_eta_j, window_eta_i) + obj.phi(C2_xi, C2_eta);
                    end
                end
            end
        end
    end
end

