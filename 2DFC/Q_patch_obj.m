classdef Q_patch_obj
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
    end
    
    methods
        function obj = Q_patch_obj(M_p, J, n_xi, n_eta, xi_start, xi_end, eta_start, eta_end, f_XY)
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
        end
        
        function phi_fh = phi(obj, xi, eta)
            xi_0 = (obj.xi_start+obj.xi_end)/2;
            eta_0 = (obj.eta_start+obj.eta_end)/2;
            R_xi = (obj.xi_start-obj.xi_end)/2;
            R_eta = (obj.eta_start-obj.eta_end)/2;
            
            phi_fh = zeros(size(xi));
            r_squared = ((xi-xi_0)/R_xi).^2 + ((eta-eta_0)/R_eta).^2;
            phi_fh(r_squared <= 1) = exp(-1./(1-r_squared(r_squared <= 1)));
        end
        
        function mesh = xi_mesh(obj)
            mesh = transpose(linspace(obj.xi_start, obj.xi_end, obj.n_xi+1));
        end
        
        function mesh = eta_mesh(obj)
            mesh = transpose(linspace(obj.eta_start, obj.eta_end, obj.n_eta+1));
        end
        
        function [XI, ETA] = xi_eta_mesh(obj)
            [XI, ETA] = meshgrid(obj.xi_mesh, obj.eta_mesh);
        end
        
        function [xi, eta, converged] = inverse_M_p(obj, x, y)
            err_guess = @(x, y, v) transpose(obj.M_p(v(1), v(2))) - [x; y];
            [v_guess, converged] = newton_solve(@(v) err_guess(x, y, v), obj.J, [0; 0], 1e-15, 100);
            xi = v_guess(1);
            eta = v_guess(2);
        end
        
        function [f_xy, phi_xy] = locally_compute(obj, x, y, d)
            % d is degree of interpolation
            [xi, eta] = obj.inverse_M_p(x, y);
            phi_xy = obj.phi(xi, eta);
            
            h_xi = (obj.xi_end-obj.xi_start)/obj.n_xi;
            h_eta = (obj.eta_end-obj.eta_start)/obj.n_eta;
            
            % j to the immediate left of point
            xi_j = floor((xi-obj.xi_start)/h_xi);
            eta_j = floor((eta-obj.eta_start)/h_eta);
            
            if (xi > obj.xi_end || xi < obj.xi_start) || (eta > obj.eta_end || eta < obj.eta_start)
                f_xy = nan;
                phi_xy = nan;
                return
            end
            
            half_d =  floor((d+1)/2);
            if mod(d, 2) ~= 0
                interpol_xi_j_mesh = transpose(xi_j-half_d+1:xi_j+half_d);
                interpol_eta_j_mesh = transpose(eta_j-half_d+1:eta_j+half_d);
            else
                interpol_xi_j_mesh = transpose(xi_j-half_d:xi_j+half_d);
                interpol_eta_j_mesh = transpose(eta_j-half_d:eta_j+half_d);
            end

            interpol_xi_j_mesh = shift_idx_mesh(interpol_xi_j_mesh, 0, obj.n_xi+1);
            interpol_eta_j_mesh = shift_idx_mesh(interpol_eta_j_mesh, 0, obj.n_eta+1);
            
            interpol_xi_mesh = h_xi*interpol_xi_j_mesh + obj.xi_start;
            interpol_eta_mesh = h_eta*interpol_eta_j_mesh + obj.eta_start;
            
            interpol_eta_exact = zeros(d+1, 1);
            for vert_idx = 1:d+1
                interpol_eta_exact(vert_idx) = barylag([interpol_eta_mesh,  obj.f_XY(interpol_eta_j_mesh+1, interpol_xi_j_mesh(vert_idx)+1)], eta);
            end
            
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

