classdef IE_curve_obj < handle
    %IE_CURVE_OBJ Curve obj for integral equation solver
    
    properties
        l_1
        l_2
        l_1_prime
        l_2_prime
        l_1_dprime
        l_2_dprime
        
        w
        w_prime
        
        curve_idx
        next_curve
    end
    
    methods
        function obj = IE_curve_obj(curve, curve_idx, p)
            %IE_CURVE_OBJ Construct an instance of this class
            %   Detailed explanation goes here
            obj.l_1 = curve.l_1;
            obj.l_2 = curve.l_2;
            obj.l_1_prime = curve.l_1_prime;
            obj.l_2_prime = curve.l_2_prime;
            obj.l_1_dprime = curve.l_1_dprime;
            obj.l_2_dprime = curve.l_2_dprime;
                        
            obj.curve_idx = curve_idx;
            
            v = @(s) (1/p-1/2)*(1-2*s).^3+1/p*(2*s-1)+1/2;
            v_prime = @(s) 2/p-6*(1/p-1/2)*(1-2*s).^2;

            obj.w = @(s) v(s).^p./(v(s).^p+v(1-s).^p);
            obj.w_prime = @(s) p*((v_prime(s).*v(s).^(p-1))./(v(s).^p+v(1-s).^p)-(v(s).^(p-1).*v_prime(s)-v(1-s).^(p-1).*v_prime(1-s)).*v(s).^p./(v(s).^p+v(1-s).^p).^2);
        end
        
        function [s_mesh, ds] = s_mesh(obj, n)
            ds = 1/n;
            s_mesh = linspace(0, 1-ds, n)';
        end
        
        function theta_mesh = theta_mesh(obj, n)
            [s_mesh, ~] = obj.s_mesh(n);
            theta_mesh = obj.w(s_mesh);
        end
        
        function u_num_curve = u_num_curve(obj, x, y, curve_param, gr_phi)
            curve_n = curve_param.curve_n(obj.curve_idx);
            start_idx = curve_param.start_idx(obj.curve_idx);
            end_idx = curve_param.end_idx(obj.curve_idx);
            [s_mesh, ds] =  obj.s_mesh(curve_n);
            theta_mesh = obj.w(s_mesh);

            u_num_curve = transpose(obj.K_general([x; y], theta_mesh).*gr_phi(start_idx:end_idx).*sqrt(obj.l_1_prime(theta_mesh).^2+obj.l_2_prime(theta_mesh).^2)) * ones(curve_n, 1) * ds;
        end
        
        function K_general_eval = K_general(obj, x, theta)
            K_general_eval = -1 / (2 * pi) * ( ...
                (x(1) - obj.l_1(theta)) .* obj.l_2_prime(theta) - ...
                (x(2) - obj.l_2(theta)) .* obj.l_1_prime(theta) ...
            ) ./ ( ...
                sqrt(obj.l_1_prime(theta).^2 + obj.l_2_prime(theta).^2) .* ...
                ( ...
                    (x(1) - obj.l_1(theta)).^2 + ...
                    (x(2) - obj.l_2(theta)).^2 ...
                ) ...
            );
        end
        
        function  K_boundary_same_point_eval = K_boundary_same_point(obj, theta)
            K_boundary_same_point_eval = -1 / (4 * pi) * ( ...
                obj.l_2_prime(theta) .* obj.l_1_dprime(theta) - ...
                obj.l_1_prime(theta) .* obj.l_2_dprime(theta) ...
            ) ./ ( ...
                (obj.l_1_prime(theta).^2 + obj.l_2_prime(theta).^2).^(3/2) ...
            );
        end
        
        function K_boundary_eval = K_boundary(obj, theta_1, theta_2, curve_2)
            K_boundary_eval = -1 / (2 * pi) * ( ...
                (obj.l_1(theta_1) - curve_2.l_1(theta_2)) .* curve_2.l_2_prime(theta_2) - ...
                (obj.l_2(theta_1) - curve_2.l_2(theta_2)) .* curve_2.l_1_prime(theta_2) ...
            ) ./ ( ...
                sqrt(curve_2.l_1_prime(theta_2).^2 + curve_2.l_2_prime(theta_2).^2) .* ...
                ( ...
                    (obj.l_1(theta_1) - curve_2.l_1(theta_2)).^2 + ...
                    (obj.l_2(theta_1) - curve_2.l_2(theta_2)).^2 ...
                ) ...
            );
        end
    end

end


