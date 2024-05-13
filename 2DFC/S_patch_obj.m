classdef S_patch_obj < Q_patch_obj
    %S_PATCH_OBJ Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        M_p_general % function handle of M_p with H as a parameter
        J_general
        h
    end
    
    methods
        function obj = S_patch_obj(M_p_general, J_general, n_xi, n_eta, xi_start, xi_end, eta_start, eta_end, f_XY, h)
            %S_PATCH_OBJ Construct an instance of this class
            %   Detailed explanation goes here
            h_eta = (eta_end-eta_start)/n_eta;
            H = h / h_eta;
            
            obj = obj@Q_patch_obj(@(xi, eta) M_p_general(xi, eta, H), @(v) J_general(v, H), n_xi, n_eta, xi_start, xi_end, eta_start, eta_end, f_XY);
            obj.M_p_general = M_p_general;
            obj.J_general = J_general;
            obj.h = h;
        end
        
        function phi_fh = phi(obj, xi, eta)
            xi_0 = (obj.xi_start+obj.xi_end)/2;
            R_xi = (obj.xi_end-obj.xi_start);

            in_rectangle = xi <= obj.xi_end & xi >= obj.xi_start & eta <= obj.eta_end & eta >= obj.eta_start;
            phi_fh = zeros(size(xi));
            
            phi_fh(in_rectangle) = exp(-1./(1-(4/R_xi.^2).*(xi(in_rectangle)-xi_0).^2));
        end
        
        function phi_fh = window_phi(obj, xi, eta)
            % one sided phi function
            xi_0 = obj.xi_end;
            R_xi = 2*(obj.xi_end-obj.xi_start);

            in_rectangle = xi <= obj.xi_end & xi >= obj.xi_start & eta <= obj.eta_end & eta >= obj.eta_start;
            phi_fh = zeros(size(xi));
            
            phi_fh(in_rectangle) = exp(-1./(1-(4/R_xi.^2).*(xi(in_rectangle)-xi_0).^2));
        end
        
        function S_fcont_patch = FC(obj, C, d, A, Q, phi_normalization)
            h_eta = (obj.eta_end-obj.eta_start)/obj.n_eta;
            
            [XI, ETA] = obj.xi_eta_mesh();
            if ~isnan(phi_normalization)
                fcont = fcont_gram_blend_S(obj.f_XY.*obj.window_phi(XI, ETA)./phi_normalization, d, A, Q);
            else
                fcont = fcont_gram_blend_S(obj.f_XY, d, A, Q);
            end
            
            S_fcont_patch = S_patch_obj(obj.M_p_general, obj.J_general, obj.n_xi, C+obj.n_eta, obj.xi_start, obj.xi_end, obj.eta_start-C*h_eta, obj.eta_end, fcont, obj.h);
        end
    end
end

