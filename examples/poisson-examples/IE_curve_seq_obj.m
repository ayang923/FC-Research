classdef IE_curve_seq_obj < handle
    %IE_CURVE_SEQ_OBJ Obj for integral equation solver
    
    properties
        n_curves
        first_curve
        last_curve
    end
    
    methods
        function obj = IE_curve_seq_obj(curve_seq, p)
            %IE_CURVE_SEQ_OBJ Construct an instance of this class
            
            obj.n_curves = curve_seq.n_curves;
            obj.first_curve = IE_curve_obj(curve_seq.first_curve, 1, p);
            
            curr_curve = curve_seq.first_curve;
            curr_ie_curve = obj.first_curve;
            for i = 2:obj.n_curves
                curr_curve= curr_curve.next_curve;
                curr_ie_curve.next_curve = IE_curve_obj(curr_curve, i, p);
                curr_ie_curve = curr_ie_curve.next_curve;
            end
            obj.last_curve = curr_ie_curve;
            curr_ie_curve.next_curve = obj.first_curve;
        end
        
        function [A, b] = construct_A_b(obj, curve_param, f)
            curve_n = curve_param.curve_n;
            n_total = curve_param.n_total;
            start_idx = curve_param. start_idx;
            end_idx = curve_param.end_idx;
            
            A = zeros(n_total, n_total);
            b = zeros(n_total, 1);


            curr_y = obj.first_curve;
            for i = 1:obj.n_curves
                [s_mesh_y, ds_y] = curr_y.s_mesh(curve_n(i));
                theta_mesh_y = curr_y.w(s_mesh_y);
                b(start_idx(i):end_idx(i)) = curr_y.w_prime(s_mesh_y).*f(curr_y.l_1(theta_mesh_y), curr_y.l_2(theta_mesh_y));

                curr_x = obj.first_curve;
                for j = 1:obj.n_curves
                    [s_mesh_x, ds_x] = curr_x.s_mesh(curve_n(j));
                    theta_mesh_x = curr_x.w(s_mesh_x);
                    [theta_2, theta_1] = meshgrid(theta_mesh_y, theta_mesh_x);

                    A_local = curr_x.w_prime(s_mesh_x).*curr_x.K_boundary(theta_1, theta_2, curr_y).*sqrt(curr_y.l_1_prime(theta_mesh_y').^2+curr_y.l_2_prime(theta_mesh_y').^2)*ds_y;

                    if i == j
                        msk_diagonals = logical(diag(ones(curve_n(j), 1)));
                        A_local(msk_diagonals) = curr_x.w_prime(s_mesh_x).*curr_x.K_boundary_same_point(theta_mesh_y).*sqrt(curr_x.l_1_prime(theta_mesh_x).^2+curr_x.l_2_prime(theta_mesh_x).^2)*ds_x + 1/2;
                    end

                    A(start_idx(j):end_idx(j), start_idx(i):end_idx(i)) = A_local;

                    curr_x = curr_x.next_curve;
                end

                curr_y = curr_y.next_curve;
            end
        end

        function u_num = u_num(obj, x, y, curve_param, gr_phi)
           curr = obj.first_curve;
           u_num = 0;
           for i = 1:obj.n_curves
                u_num =  u_num + curr.u_num_curve(x, y, curve_param, gr_phi);
                curr = curr.next_curve;
            end
        end
    end
end

