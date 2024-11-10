classdef Curve_seq_obj
    properties
        first_curve
        n_curves
    end
    
    methods
        function obj = Curve_seq_obj(first_curve, n_curves)
            obj.first_curve = first_curve;
            obj.n_curves = n_curves;
        end
        
        function patches = construct_patches(obj, f, d, eps_xi_eta, eps_xy)
            curr = obj.first_curve;
            prev_S_patch = nan;
            prev_C_patch = nan;
            patches = cell(obj.n_curves*2, 1);
            
            figure;
            for i = 1:obj.n_curves
                S_patch = curr.construct_S_patch(f, d, eps_xi_eta, eps_xy);
                C_patch = curr.construct_C_patch(f, d, eps_xi_eta, eps_xy);
                
                [X, Y] = S_patch.Q.xy_mesh;
                scatter(X(:), Y(:));%,patches{2*i-1}.Q.f_XY(:));
                hold on;
                [X, Y] = C_patch.L.xy_mesh;
                scatter(X(:), Y(:));%,patches{2*i}.L.f_XY);
                [X, Y] = C_patch.W.xy_mesh;
                scatter(X(:), Y(:));
                

                if i ~= 1
                    prev_C_patch.apply_w_W(prev_S_patch)
                    prev_C_patch.apply_w_L(S_patch)
                end
                
                patches{2*i-1} = S_patch;
                patches{2*i} = C_patch;
                
                prev_S_patch = S_patch;
                prev_C_patch = C_patch;
                curr = curr.next_curve;
            end
            C_patch.apply_w_W(S_patch);
            C_patch.apply_w_L(patches{1});
        end
        
        function [boundary_X, boundary_Y] = construct_boundary_mesh(obj, n_r)
            curr = obj.first_curve;
            n_points = 0;
            for i = 1:obj.n_curves
                n_points = n_points + (curr.n-1)*n_r + 1;
                curr = curr.next_curve;
            end
            
            boundary_X = zeros(n_points, 1);
            boundary_Y = zeros(n_points, 1);
            curr_idx = 1;
            curr = obj.first_curve;
            for i = 1:obj.n_curves
                boundary_X(curr_idx:curr_idx+(curr.n-1)*n_r) = curr.l_1(linspace(0, 1, (curr.n-1)*n_r + 1)');
                boundary_Y(curr_idx:curr_idx+(curr.n-1)*n_r) = curr.l_2(linspace(0, 1, (curr.n-1)*n_r + 1)');
                
                curr = curr.next_curve;
            end
        end
    end
end

