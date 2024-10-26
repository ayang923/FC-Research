classdef C2_patch_obj < handle
    %C2_PATCH_OBJ Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        W
        L
    end
    
    methods
        function obj = C2_patch_obj(M_p, J, eps_xi_eta, eps_xy, n_xi, n_eta, d, f_W, f_L)
            %C2_PATCH_OBJ Construct an instance of this class
            %   Detailed explanation goes here
            h_xi = 1./(n_xi-1);
            h_eta = 1./(n_eta-1);
            obj.W = Q_patch_obj(M_p, J, eps_xi_eta, eps_xy, n_xi, d, 0, 1, 0, (d-1)*h_eta, f_W, nan);            
            obj.L = Q_patch_obj(M_p, J, eps_xi_eta, eps_xy, d, n_eta-d+1, 0, (d-1)*h_xi, (d-1)*h_eta, 1, f_L, nan);
        end
                
        function [C2_fcont_patch_xi, C2_fcont_patch_eta, C2_fcont_patch_corner] = FC(obj, C, n_r, d, A, Q, phi_normalization)
            [h_xi, h_eta] = obj.h_mesh;
            [XI, ETA] = obj.xi_eta_mesh();
            
            if ~isnan(phi_normalization)
                fcont = fcont_gram_blend_C2(obj.f_XY.*obj.phi(XI, ETA)./phi_normalization, d, A, Q);
            else
                fcont = fcont_gram_blend_C2(obj.f_XY, d, A, Q);
            end
            
            C2_fcont_patch_xi =  Q_patch_obj(obj.M_p, obj.J, obj.eps_xi_eta, obj.eps_xy, obj.n_xi, C*n_r+1, obj.xi_start, obj.xi_end, obj.eta_start-(C)*h_eta, obj.eta_start, fcont(1:C*n_r+1, C*n_r+1:end), obj.phi);
            C2_fcont_patch_eta = Q_patch_obj(obj.M_p, obj.J, obj.eps_xi_eta, obj.eps_xy, C*n_r+1, obj.n_eta, obj.xi_start-(C)*h_xi, obj.xi_start, obj.eta_start, obj.eta_end, fcont(C*n_r+1:end, 1:C*n_r+1), obj.phi);
            C2_fcont_patch_corner = Q_patch_obj(obj.M_p, obj.J, obj.eps_xi_eta, obj.eps_xy, C*n_r+1, C*n_r+1, obj.xi_start-(C)*h_xi, obj.xi_start, obj.eta_start-C*h_eta, obj.eta_start, fcont(1:C*n_r+1, 1:C*n_r+1), obj.phi);
        end
        
        function [C2_W_norm, C2_L_norm, xi_norm, eta_norm] = compute_phi_normalization(obj, window_patch_xi, window_patch_eta)
            [C2_W_norm, xi_norm] = obj.W.compute_phi_normalization_xi_right(window_patch_xi.Q_patch);
            [C2_L_norm, eta_norm] = obj.L.compute_phi_normalization_eta_up(window_patch_eta.Q_patch);
        end
    end
end





