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
                
        function [C2_fcont_patch_xi, C2_fcont_patch_eta, C2_fcont_patch_corner] = FC(obj, C, n_r, d, A, Q, W_norm, L_norm)
            [h_xi, h_eta] = obj.W.h_mesh; %h mesh is the same for both W and L
            if ~isnan(W_norm)
                [XI_W, ETA_W] = obj.W.xi_eta_mesh;
                f_XY_W = obj.W.f_XY .* obj.W.phi(XI_W, ETA_W) ./ W_norm;
            else
                f_XY_W = obj.W.f_XY;
            end
            
            if ~isnan(L_norm)
                [XI_L, ETA_L] = obj.L.xi_eta_mesh;
                f_XY_L = obj.L.f_XY .* obj.L.phi(XI_L, ETA_L) ./ L_norm;
            else
                f_XY_L = obj.L.f_XY;
            end
            
            fcont_W = fcont_gram_blend_S(f_XY_W, d, A, Q);
            fcont_corner = transpose(fcont_gram_blend_S(fcont_W', d, A, Q));
            fcont_L = transpose(fcont_gram_blend_S([f_XY_W(:, 1:d); f_XY_L(2:end, :)]', d, A, Q));
            
            
            C2_fcont_patch_xi =  Q_patch_obj(obj.W.M_p, obj.W.J, obj.W.eps_xi_eta, obj.W.eps_xy, obj.W.n_xi, C*n_r+1, 0, 1, -(C)*h_eta, 0, fcont_W, nan);
            C2_fcont_patch_eta = Q_patch_obj(obj.L.M_p, obj.L.J, obj.L.eps_xi_eta, obj.L.eps_xy, C*n_r+1, obj.L.n_eta + d-1, -(C)*h_xi, 0, 0, 1, fcont_L, nan);
            C2_fcont_patch_corner = Q_patch_obj(obj.W.M_p, obj.W.J, obj.W.eps_xi_eta, obj.W.eps_xy, C*n_r+1, C*n_r+1, -(C)*h_xi, 0, -C*h_eta, 0, fcont_corner, nan);
        end
        
        function [C2_W_norm, C2_L_norm, xi_norm, eta_norm] = compute_phi_normalization(obj, window_patch_xi, window_patch_eta)
            [C2_W_norm, xi_norm] = obj.W.compute_phi_normalization_xi_right(window_patch_xi.Q_patch);
            [C2_L_norm, eta_norm] = obj.L.compute_phi_normalization_eta_up(window_patch_eta.Q_patch);
        end
    end
end





