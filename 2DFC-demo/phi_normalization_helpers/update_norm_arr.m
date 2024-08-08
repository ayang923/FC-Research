function [norm_arr] = update_norm_arr(norm_arr, window_patch, overlap_X, overlap_Y, overlap_XI_j, overlap_ETA_j, initial_guesses)
    for i = 1:size(overlap_X, 1)
        if mod(i, 2) == 1
            j_lst = 1:size(overlap_X, 2);
        else
            j_lst = size(overlap_X, 2):-1:1;
        end
        
        for j = j_lst
            [window_patch_xi, window_patch_eta, converged] = window_patch.inverse_M_p(overlap_X(i, j), overlap_Y(i, j), initial_guesses);
            
            if converged && window_patch.in_patch(window_patch_xi, window_patch_eta)
                xi_i = overlap_XI_j(i, j) + 1;
                xi_j = overlap_ETA_j(i, j) + 1;
                norm_arr(xi_j, xi_i) = norm_arr(xi_j, xi_i) + window_patch.phi(window_patch_xi, window_patch_eta);
            elseif ~converged
                warning("Nonconvergence in computing C2_norm")
            end
            
            if converged
                initial_guesses = [window_patch_xi; window_patch_eta]; % using previous point as next initial guess
            end
        end
    end
end
