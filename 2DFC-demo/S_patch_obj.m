classdef S_patch_obj < handle
    %S_PATCH_OBJ Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Q_patch
        M_p_general % function handle of M_p with H as a parameter
        J_general
        h
    end
    
    methods
        function obj = S_patch_obj(M_p_general, J_general, h, eps_xi_eta, eps_xy, n_xi, n_eta, f_XY)
            %S_PATCH_OBJ Construct an instance of this class
            %   Detailed explanation goes here
            h_eta = 1/(n_eta-1);
            H = h / h_eta;
            
            obj.Q_patch = Q_patch_obj(@(xi, eta) M_p_general(xi, eta, H), @(v) J_general(v, H), eps_xi_eta, eps_xy, n_xi, n_eta, 0, 1, 0, 1, f_XY, nan);
            
            obj.M_p_general = M_p_general;
            obj.J_general = J_general;
            obj.h = h;
        end
        
        function S_fcont_patch = FC(obj, C, n_r, d, A, Q, phi_normalization)
            h_eta = (obj.Q_patch.eta_end-obj.Q_patch.eta_start)/(obj.Q_patch.n_eta-1);
            
            [XI, ETA] = obj.Q_patch.xi_eta_mesh();
            if ~isnan(phi_normalization)
                fcont = fcont_gram_blend_S(obj.Q_patch.f_XY.*obj.Q_patch.phi(XI, ETA)./phi_normalization, d, A, Q);
            else
                fcont = fcont_gram_blend_S(obj.Q_patch.f_XY, d, A, Q);
            end
            
            S_fcont_patch = Q_patch_obj(obj.Q_patch.M_p, obj.Q_patch.J, obj.Q_patch.eps_xi_eta, obj.Q_patch.eps_xy, obj.Q_patch.n_xi, C*n_r+1, obj.Q_patch.xi_start, obj.Q_patch.xi_end, obj.Q_patch.eta_start-C*h_eta, obj.Q_patch.eta_start, fcont, nan);
        end
    end
end

