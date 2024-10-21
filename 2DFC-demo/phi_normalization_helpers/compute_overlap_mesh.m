function [XI_overlap, ETA_overlap, XI_j, ETA_j] = compute_overlap_mesh(main_patch, xi_eta_corner, quadrant_overlap)
    [h_xi, h_eta] = main_patch.h_mesh();
    j_corner =  (xi_eta_corner - [main_patch.xi_start; main_patch.eta_start])./ [h_xi; h_eta];
    
    if quadrant_overlap == 1
        j_corner_grid = [ceil(j_corner(1)); ceil(j_corner(2))];
        [XI_j, ETA_j] = meshgrid(max([j_corner_grid(1); 0]):main_patch.n_xi-1, max([j_corner_grid(2); 0]):main_patch.n_eta-1);
    elseif quadrant_overlap == 2
        j_corner_grid = [floor(j_corner(1)); ceil(j_corner(2))];
        [XI_j, ETA_j] = meshgrid(0:min([main_patch.n_xi-1; j_corner_grid(1)]), max([j_corner_grid(2); 0]):main_patch.n_eta-1);
    elseif quadrant_overlap == 3
        j_corner_grid = [floor(j_corner(1)); floor(j_corner(2))];
        [XI_j, ETA_j] = meshgrid(0:min([main_patch.n_xi-1; j_corner_grid(1)]), 0:min([main_patch.n_eta-1; j_corner_grid(2)]));
    elseif quadrant_overlap == 4
        j_corner_grid = [ceil(j_corner(1)); floor(j_corner(2))];
        [XI_j, ETA_j] = meshgrid(max([j_corner_grid(1); 0]):main_patch.n_xi-1, 0:min([main_patch.n_eta-1; j_corner_grid(2)]));
    else
        error("Invalid Quadrant Number")
    end
    XI_overlap = XI_j * h_xi + main_patch.xi_start;
    ETA_overlap = ETA_j * h_eta + main_patch.eta_start;
end

