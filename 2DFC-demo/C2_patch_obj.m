classdef C2_patch_obj < Q_patch_obj
    %C2_PATCH_OBJ Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
        function obj = C2_patch_obj(M_p, J, eps_xi_eta, eps_xy, n_xi, n_eta, xi_start, xi_end, eta_start, eta_end, f_XY, phi)
            %C2_PATCH_OBJ Construct an instance of this class
            %   Detailed explanation goes here
            obj = obj@Q_patch_obj(M_p, J, eps_xi_eta, eps_xy, n_xi, n_eta, xi_start, xi_end, eta_start, eta_end, f_XY, phi);
        end
                
        function [C2_fcont_patch_xi, C2_fcont_patch_eta, C2_fcont_patch_corner] = FC(obj, C, n_r, d, A, Q, phi_normalization)
            h_xi = (obj.xi_end-obj.xi_start)/obj.n_xi;
            h_eta = (obj.eta_end-obj.eta_start)/obj.n_eta;
            
            [XI, ETA] = obj.xi_eta_mesh();
            
            if ~isnan(phi_normalization)
                fcont = fcont_gram_blend_C2(obj.f_XY.*obj.phi(XI, ETA)./phi_normalization, d, A, Q);
            else
                fcont = fcont_gram_blend_C2(obj.f_XY, d, A, Q);
            end
            
            C2_fcont_patch_xi =  C2_patch_obj(obj.M_p, obj.J, obj.eps_xi_eta, obj.eps_xy, obj.n_xi, C*n_r, obj.xi_start, obj.xi_end, obj.eta_start-(C)*h_eta, obj.eta_start, fcont(1:C*n_r+1, C*n_r+1:end), obj.phi);
            C2_fcont_patch_eta = C2_patch_obj(obj.M_p, obj.J, obj.eps_xi_eta, obj.eps_xy, C*n_r, obj.n_eta, obj.xi_start-(C)*h_xi, obj.xi_start, obj.eta_start, obj.eta_end, fcont(C*n_r+1:end, 1:C*n_r+1), obj.phi);
            C2_fcont_patch_corner = C2_patch_obj(obj.M_p, obj.J, obj.eps_xi_eta, obj.eps_xy, C*n_r, C*n_r, obj.xi_start-(C)*h_xi, obj.xi_start, obj.eta_start-C*h_eta, obj.eta_start, fcont(1:C*n_r+1, 1:C*n_r+1), obj.phi);
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
            
            window_patch_xi.phi = @(xi, eta) window_patch_xi.phi_1D((xi-window_xi_xi_0)/window_xi_R_xi);
            
            window_eta_R_xi = window_eta_corner(1) - window_patch_eta.xi_start;
            window_eta_xi_0 = window_patch_eta.xi_start;
            
            window_patch_eta.phi = @(xi, eta) window_patch_eta.phi_1D((xi-window_eta_xi_0)/window_eta_R_xi) ;
            
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
            
            xi_norm = update_norm_arr(xi_norm, obj, window_xi_overlap_X, window_xi_overlap_Y, window_xi_overlap_XI_j, window_xi_overlap_ETA_j, nan);
            
            [eta_window_XI, eta_window_ETA] = window_patch_eta.xi_eta_mesh();
            eta_norm = window_patch_eta.phi(eta_window_XI, eta_window_ETA);
                        
            [window_eta_overlap_XI, window_eta_overlap_ETA, window_eta_overlap_XI_j, window_eta_overlap_ETA_j] = compute_overlap_mesh(window_patch_eta, window_eta_corner, 3);
            [window_eta_overlap_X, window_eta_overlap_Y] = window_patch_eta.convert_to_XY(window_eta_overlap_XI, window_eta_overlap_ETA);
            
            eta_norm = update_norm_arr(eta_norm, obj, window_eta_overlap_X, window_eta_overlap_Y, window_eta_overlap_XI_j, window_eta_overlap_ETA_j, nan);
        end
    end
end





