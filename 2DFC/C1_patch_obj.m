classdef C1_patch_obj < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        B
        B_c % B complement
        n_bound
        phi
        h
    end
    
    methods
        function obj = C1_patch_obj(M_p, J, n_bound, phi, f_B, f_B_c)
            obj.B = Q_patch_obj(M_p, J, n_bound, n_bound*2, 1/2, 1, 0, 1, f_B, phi);
            obj.B_c = Q_patch_obj(M_p, J, n_bound, n_bound, 0, 1/2, 1/2, 1, f_B_c, phi);
            
            obj.h = 1/(2*n_bound);
            
            obj.n_bound = n_bound;
            obj.phi = phi;
        end
        
        function [boundary_mesh_xi, boundary_mesh_eta] = boundary_mesh(obj)
            n_boundary = 10;
            boundary_mesh_xi = [zeros(n_boundary, 1); linspace(0, 1, n_boundary)'; ones(n_boundary, 1); linspace(1, 1/2, n_boundary)'; 1/2*ones(n_boundary, 1); linspace(1/2, 0, n_boundary)'];
            boundary_mesh_eta = [linspace(1/2, 1, n_boundary)'; ones(n_boundary, 1); linspace(1, 0, n_boundary)'; zeros(n_boundary, 1); linspace(0, 1/2, n_boundary)'; 1/2*ones(n_boundary, 1)];
        end
        
        function [boundary_mesh_x, boundary_mesh_y] = boundary_mesh_xy(obj)
            [boundary_mesh_xi, boundary_mesh_eta] = obj.boundary_mesh();
            [boundary_mesh_x, boundary_mesh_y] = obj.B.convert_to_XY(boundary_mesh_xi, boundary_mesh_eta);
        end
        
        function  [C1_fcont_patch_xi, C1_fcont_patch_eta_refined, C1_fcont_patch_eta_unrefined, B_fcont, B_c_minus_fcont] = FC(obj, C, n_r, d, A, Q, f, M)
            B_fcont = transpose(fcont_gram_blend_S(obj.B.f_XY', d, A, Q));
            B_fc_unrefined = B_fcont(:, 1:n_r:n_r*C+1); % includes boundary point from where continuation is done
            
            C1_fcont_patch_xi = Q_patch_obj(obj.B.M_p, obj.B.J, C*n_r, obj.n_bound, obj.B.xi_start-C*obj.h, obj.B.xi_start, obj.B.eta_start, obj.B_c.eta_start, B_fcont(1:obj.n_bound+1, 1:C*n_r+1), obj.B.phi);
            
            if obj.n_bound < C
                B_c_minus_fcont = padarray(obj.B_c.f_XY, [0 C - obj.n_bound], 0, 'pre') - B_fc_unrefined(obj.n_bound+1:end, :);
                B_c_minus_fcont_xi_start = obj.B_c.xi_start - (C-obj.n_bound)*obj.h;
            else
                [B_c_unrefined_f_XY, B_c_refined_f_XY] = obj.refine_B_c(C, n_r, f, M);
%                 B_c_minus_fcont = obj.B_c.f_XY - padarray(B_fc_unrefined(obj.n_bound+1:end, :), [0 obj.n_bound-C], 0, 'pre');
                B_c_minus_fcont = B_c_refined_f_XY - B_fcont(obj.n_bound+1:end, 1:n_r*C+1);
                
                B_c_fcont_refined = fcont_gram_blend_S(B_c_minus_fcont, d, A, Q);
                B_c_fcont_unrefined = fcont_gram_blend_S(B_c_unrefined_f_XY, d, A, Q);
                
                C1_fcont_patch_eta_refined = Q_patch_obj(obj.B_c.M_p, obj.B_c.J, C*n_r, C*n_r, obj.B_c.xi_end-C*obj.h, obj.B_c.xi_end, obj.B_c.eta_start-C*obj.h, obj.B_c.eta_start, B_c_fcont_refined(1:C*n_r+1, :), obj.B_c.phi);
                C1_fcont_patch_eta_unrefined = Q_patch_obj(obj.B_c.M_p, obj.B_c.J, obj.n_bound-C, C*n_r, obj.B_c.xi_start, obj.B_c.xi_end-C*obj.h, obj.B_c.eta_start-C*obj.h, obj.B_c.eta_start, B_c_fcont_unrefined(1:C*n_r+1, :), obj.B_c.phi);
            end

%             B_c_fcont = fcont_gram_blend_S(B_c_minus_fcont, d, A, Q);
%             C1_fcont_patch_eta = Q_patch_obj(obj.B_c.M_p, obj.B_c.J, max(C, obj.n_bound), C*n_r, B_c_minus_fcont_xi_start, obj.B_c.xi_end, obj.B_c.eta_start-C*obj.h, obj.B_c.eta_start, B_c_fcont(1:C*n_r+1, :), obj.B_c.phi);
            
        end
        
        function [B_c_unrefined_f_XY, B_c_refined_f_XY] = refine_B_c(obj, C, n_r, f, M)
            B_c_unrefined_f_XY = obj.B_c.f_XY(:, 1:obj.n_bound-C+1);
            
            B_c_refined_xi_mesh = transpose((obj.B_c.xi_end - C*obj.h):(obj.h/n_r):obj.B_c.xi_end);
                                    
            B_c_refined_f_XY = zeros(obj.B_c.n_eta+1, length(B_c_refined_xi_mesh));
            half_M = floor((M+1)/2);
            for eta_j = 0:obj.n_bound
                for xi_j = obj.n_bound-C:(obj.n_bound-1)
                    if mod(M, 2) ~= 0
                        interpol_xi_j_mesh = transpose(xi_j-half_M+1:xi_j+half_M);
                    else
                        interpol_xi_j_mesh = transpose(xi_j-half_M:xi_j+half_M);
                    end

                    interpol_xi_j_mesh = shift_idx_mesh(interpol_xi_j_mesh, 0, obj.n_bound);
                    interpol_xi_mesh = obj.h*interpol_xi_j_mesh + obj.B_c.xi_start;
                    interpol_val = obj.B_c.f_XY(eta_j+1, interpol_xi_j_mesh+1)';
                    
                    xi_eval_j = ((xi_j-(obj.n_bound-C))*n_r):((xi_j-(obj.n_bound-C)+1)*n_r);
                    
                    B_c_refined_f_XY(eta_j+1, xi_eval_j+1) = barylag([interpol_xi_mesh, interpol_val], B_c_refined_xi_mesh(xi_eval_j+1));
                end
            end
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

