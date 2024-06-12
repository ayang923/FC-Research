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
        
%         function phi_fh = phi(obj, xi, eta)
%             xi_0 = (obj.xi_start+obj.xi_end)/2;
%             eta_0 = (obj.eta_start+obj.eta_end)/2;
%             R_xi = (obj.xi_start-obj.xi_end);
%             R_eta = (obj.eta_start-obj.eta_end);
% 
%             in_rectangle = xi <= obj.xi_end & xi >= obj.xi_start & eta <= obj.eta_end & eta >= obj.eta_start;
%             phi_fh = zeros(size(xi));
%             
%             phi_fh(in_rectangle) = exp(-1./(1-(4/R_xi.^2).*(xi(in_rectangle)-xi_0).^2)).*exp(-1./(1-(4/R_eta.^2).*(eta(in_rectangle)-eta_0).^2));
%         end
        
        function [h_xi, h_eta] = h_mesh(obj)
            h_xi =(obj.xi_end-obj.xi_start)./obj.n_xi;
            h_eta = (obj.eta_end-obj.eta_start) ./ obj.n_eta;
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
        
        function [X, Y] = xy_mesh(obj)
            [XI, ETA] = obj.xi_eta_mesh();
            xy = obj.M_p(XI(:), ETA(:));
            
            X = reshape(xy(:, 1), size(XI));
            Y = reshape(xy(:, 2), size(XI));
        end
        
        function patch_msk = in_patch(obj, xi, eta)
            eps = 1e-15;
            patch_msk = xi >= obj.xi_start-eps & xi <= obj.xi_end+eps & eta >= obj.eta_start-eps & eta <= obj.eta_end+eps;
        end
        
        function patch_msk = in_patch_exact(obj, xi, eta)
            patch_msk =  xi >= obj.xi_start & xi <= obj.xi_end & eta >= obj.eta_start & eta <= obj.eta_end;
        end
        
        function [xi, eta, converged] = inverse_M_p(obj, x, y)
            % need to change randomization, not giving  consistently good
            % results, probably ocnvergence issues.
            N = 15;
            xi_initial = unifrnd(obj.xi_start, obj.xi_end, N, 1);
            eta_initial = unifrnd(obj.eta_start, obj.eta_end, N, 1);
            
            err_guess = @(x, y, v) transpose(obj.M_p(v(1), v(2))) - [x; y];
            eps =  1e-15;
            for k = 1:N
                initial_guess = [xi_initial(k); eta_initial(k)];
                [v_guess, converged] = newton_solve(@(v) err_guess(x, y, v), obj.J, initial_guess,eps, 100);
                
                xi = v_guess(1);
                eta = v_guess(2);
                if converged && obj.in_patch(xi, eta)
                    return
                end
            end
        end
        
        function [f_xy, phi_xy, in_range] = locally_compute(obj, x, y, d)
            % d is degree of interpolation
            [xi, eta, converged] = obj.inverse_M_p(x, y);
            
            % check if xi, eta are within bounds of Q
            if ~converged || ~obj.in_patch(xi, eta)
                f_xy = converged;
                phi_xy = nan;
                in_range = false;
                return
            end
            
            % partition of unity value
            phi_xy = obj.phi(xi, eta);
            
            h_xi = (obj.xi_end-obj.xi_start)/obj.n_xi;
            h_eta = (obj.eta_end-obj.eta_start)/obj.n_eta;
            
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

            interpol_xi_j_mesh = shift_idx_mesh(interpol_xi_j_mesh, 0, obj.n_xi);
            interpol_eta_j_mesh = shift_idx_mesh(interpol_eta_j_mesh, 0, obj.n_eta);
            
            interpol_xi_mesh = h_xi*interpol_xi_j_mesh + obj.xi_start;
            interpol_eta_mesh = h_eta*interpol_eta_j_mesh + obj.eta_start;
            
            % first 1D interpolation
            interpol_eta_exact = zeros(d+1, 1);
            for vert_idx = 1:d+1
                interpol_eta_exact(vert_idx) = barylag([interpol_eta_mesh,  obj.f_XY(interpol_eta_j_mesh+1, interpol_xi_j_mesh(vert_idx)+1)], eta);
            end
            
            % second 1D interpolation
            f_xy = barylag([interpol_xi_mesh, interpol_eta_exact], xi);
            in_range = true;
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

