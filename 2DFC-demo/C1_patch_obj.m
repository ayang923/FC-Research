classdef C1_patch_obj < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        L
        W % B complement
    end
    
    methods
        function obj = C1_patch_obj(M_p, J, eps_xi_eta, eps_xy, n_xi, n_eta, d, f_L, f_W)
            assert(mod(n_xi, 2) == 1 && mod(n_eta, 2) == 1, "n_xi and n_eta must be odd");
            
            h_xi = 1./(n_xi-1);
            h_eta = 1./(n_eta-1);
            
            obj.L = Q_patch_obj(M_p, J, eps_xi_eta, eps_xy, d, (n_eta+1)/2+(d-1), 1/2, 1/2+(d-1)*h_xi, 0, 1/2+(d-1)*h_eta, f_L, nan);
            obj.W = Q_patch_obj(M_p, J, eps_xi_eta, eps_xy, (n_xi+1)/2, d, 0, 1/2, 1/2, 1/2+(d-1)*h_eta, f_W, nan);
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
        
        function  [C1_fcont_patch_L, C1_fcont_patch_W_refined, C1_fcont_patch_W_unrefined] = FC(obj, C, n_r, d, A, Q, M, phi_L_normalization, phi_W_normalization)
            if ~isnan(phi_L_normalization)
                [XI, ETA] = obj.L.xi_eta_mesh;
                L_f_XY = obj.L.f_XY .* obj.L.phi(XI, ETA) ./ phi_L_normalization;
            else
                L_f_XY = obj.L.f_XY;
            end
            
            if ~isnan(phi_W_normalization)
                [XI, ETA] = obj.W.xi_eta_mesh;
                W_f_XY = obj.W.f_XY .* obj.W.phi(XI, ETA) ./ phi_W_normalization;
            else
                W_f_XY = obj.W.f_XY;
            end
            [h_xi, h_eta] = obj.L.h_mesh;

            L_fcont = transpose(fcont_gram_blend_S(L_f_XY', d, A, Q));
            C1_fcont_patch_L = Q_patch_obj(obj.L.M_p, obj.L.J, obj.L.eps_xi_eta, obj.L.eps_xy, C*n_r+1,  obj.L.n_eta - (d-1), 1/2-C*h_xi, 1/2, 0, 1/2, L_fcont(1:(obj.L.n_eta - (d-1)), :), nan);
            
            if 1/2-C*h_xi < 0
            
            else
                [W_unrefined_f_XY, W_refined_f_XY] = obj.refine_W(W_f_XY, C, n_r, M);
                W_minus_fcont = W_refined_f_XY - L_fcont(end-(d-1):end, :);
                
                W_fcont_refined = fcont_gram_blend_S(W_minus_fcont, d, A, Q);
                W_fcont_unrefined = fcont_gram_blend_S(W_unrefined_f_XY, d, A, Q);
                
                C1_fcont_patch_W_refined = Q_patch_obj(obj.W.M_p, obj.W.J, obj.W.eps_xi_eta, obj.W.eps_xy, C*n_r+1, C*n_r+1, 1/2-C*h_xi, 1/2, 1/2-C*h_eta, 1/2, W_fcont_refined, nan);
                C1_fcont_patch_W_unrefined = Q_patch_obj(obj.W.M_p, obj.W.J, obj.W.eps_xi_eta, obj.W.eps_xy, obj.W.n_xi-C, C*n_r+1, 0, 1/2-C*h_xi, 1/2-C*h_eta, 1/2, W_fcont_unrefined, nan);
            end
        end
        
        function [W_unrefined_f_XY, W_refined_f_XY] = refine_W(obj, W_f_XY, C, n_r, M)
            [h_xi, ~] = obj.W.h_mesh;
            
            W_unrefined_f_XY = W_f_XY(:, 1:(obj.W.n_xi-C));
            
            W_refined_xi_mesh = transpose((1/2 - C*h_xi):(h_xi/n_r):1/2);
                                    
            W_refined_f_XY = zeros(obj.W.n_eta, length(W_refined_xi_mesh));
            half_M = floor(M/2);
            for eta_j = 0:obj.W.n_eta-1
                for xi_j = obj.W.n_xi-C-1:(obj.W.n_xi-2)
                    if mod(M, 2) ~= 0
                        interpol_xi_j_mesh = transpose(xi_j-half_M:xi_j+half_M);
                    else
                        interpol_xi_j_mesh = transpose(xi_j-half_M+1:xi_j+half_M);
                    end

                    interpol_xi_j_mesh = shift_idx_mesh(interpol_xi_j_mesh, 0, obj.W.n_xi-1);
                    interpol_xi_mesh = h_xi*interpol_xi_j_mesh;
                    interpol_val = W_f_XY(eta_j+1, interpol_xi_j_mesh+1)';
                    
                    xi_eval_j = ((xi_j-(obj.W.n_xi-C-1))*n_r):((xi_j-(obj.W.n_xi-C-1)+1)*n_r);
                    
                    W_refined_f_XY(eta_j+1, xi_eval_j+1) = barylag([interpol_xi_mesh, interpol_val], W_refined_xi_mesh(xi_eval_j+1));
                end
            end
        end
        
        function [C1_L_norm, C1_W_norm, xi_norm, eta_norm] = compute_phi_normalization(obj, window_patch_xi, window_patch_eta)
            [C1_W_norm, xi_norm] = obj.W.compute_phi_normalization_xi_left(window_patch_xi.Q_patch);
            [C1_L_norm, eta_norm] = obj.L.compute_phi_normalization_eta_down(window_patch_eta.Q_patch);
        end
    end
end


