classdef S_patch_obj < Q_patch_obj
    %S_PATCH_OBJ Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        M_p_general % function handle of M_p with H as a parameter
        J_general
        h
    end
    
    methods
        function obj = S_patch_obj(M_p_general, J_general, n_xi, n_eta, xi_start, xi_end, eta_start, eta_end, f_XY, h, phi)
            %S_PATCH_OBJ Construct an instance of this class
            %   Detailed explanation goes here
            h_eta = (eta_end-eta_start)/n_eta;
            H = h / h_eta;
            
            obj = obj@Q_patch_obj(@(xi, eta) M_p_general(xi, eta, H), @(v) J_general(v, H), n_xi, n_eta, xi_start, xi_end, eta_start, eta_end, f_XY, phi);
            obj.M_p_general = M_p_general;
            obj.J_general = J_general;
            obj.h = h;
        end
        
        function S_fcont_patch = FC(obj, C, n_r, d, A, Q, phi_normalization)
            h_eta = (obj.eta_end-obj.eta_start)/obj.n_eta;
            
            [XI, ETA] = obj.xi_eta_mesh();
            if ~isnan(phi_normalization)
                fcont = fcont_gram_blend_S(obj.f_XY.*obj.phi(XI, ETA)./phi_normalization, d, A, Q);
            else
                fcont = fcont_gram_blend_S(obj.f_XY, d, A, Q);
            end
            
            S_fcont_patch = S_patch_obj(obj.M_p_general, obj.J_general, obj.n_xi, C*n_r-1, obj.xi_start, obj.xi_end, obj.eta_start-C*h_eta, obj.eta_start, fcont(1:C*n_r, :), obj.h, obj.phi);
        end
    end
end

