% R_CARTESIAN_MESH_OBJ - Generic patch objects used by S-type, C1-type, and C2-type
% patches. Assumes a rectangular domain in parameter space.
%
% Properties:
%   M_p - parametrization that maps Q from parameter space to real space,
%   assumed to be vectorized
%   J - Jacobian of M_p - 
%   n_xi - number of points xi axis is discretized into
%   n_eta - number of points eta axis is discretized into
%   xi_start - minimum value of xi in Q
%   xi_end - maximum value of xi in Q
%   eta_start - minimum value of eta in Q
%   eta_end - maximum value of eta in Q
%   f_XY - function values associated with patch on discretized mesh
%   x_min - minimum x value of M_p(Q)
%   x_max - maximum x value of M_p(Q)
%   y_min - minimum y alue of M_p(Q)
%   y_max - maximum y value of M_p(Q)
%   w_1D - one dimensional window function used to construct partition of
%       unity
%   w - unnormalized partition of unity function for this patch
%   eps_xi_eta - error tolerance in xi-eta space
%   eps_xy - error tolerance in x-y space, should be dependent on
%       eps_xi_eta and the maximum value of the J in Q
%
% Methods:
%   Q_patch_obj - Class constructor.
%   h_mesh - returns meshsize for xi and eta for object
%   xi_mesh - returns discretized xi mesh for object
%   eta_mesh - returns discretized eta mesh for object
%   xi_eta_mesh - returns entire discretzied (xi, eta) mesh associated with
%       patch -- i.e. the domain of the patch
%   xy_mesh - returns M_p(xi_eta_mesh)
%   boundary_mesh - returns discretized mesh of boundary in parameter space
%   boundary_mesh_xy - returns M_p(boundary_mesh_xy)
%   convert_to_XY - converts M_p(xi, eta) for given vector/matrices M_p,
%   in_patch - returns whether a given mesh in xi-eta space is in the domain
%   round_boundary_points - rounds points within the prescribed error tolerance near the boundary to exact
%       boundary points
%   inverse_M_p - computes the inverse of M_p numerically using Newton's
%       method
%   locally_compute - computes the function value of some xi-eta point not
%       in the domain of the patch using polynomial interpolation
%   compute_w_normalization_xi_right - computes partition of unity
%       normalization values for a "main" Q patch and "window" Q patch where
%       the window Q patch is to the right of the main patch
%   compute_w_normalization_xi_left - computes partition of unity
%       normalization values where window patch is to the left of the
%       main patch
%   compute_w_normalization_eta_up - computes partition of unity
%       normalization values where window patch is above the main patch
%   compute_w_normalization_eta_down - computes partition of unity
%       normalization values where window patch is below the main patch
%
%
%   Note the compute_w_normalization functions operate on this
%   assumptions about the window patches:
%       - the window patch is parametrized such that the "bounding" edge of
%       the window patch (edge of window patch that is"contained" within
%       main patch) corresponds to eta-axis in the window patch's parameter
%       space
%
% Author: Allen Yang
% Email: aryang@caltech.edu

classdef R_cartesian_mesh_obj < handle
    %Cartesian mesh being interpolated onto
    
    properties
        x_start
        x_end
        y_start
        y_end
        h
        n_x
        n_y
        x_mesh
        y_mesh
        R_X
        R_Y
        R_idxs
        
        boundary_X
        boundary_Y
        in_interior
        interior_idxs
        f_R
        
        fc_coeffs
    end
    
    methods
        function obj = R_cartesian_mesh_obj(x_start,x_end, y_start, y_end, h, boundary_X, boundary_Y)
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here
            obj.x_start = x_start;
            obj.y_start = y_start;
            obj.x_end = ceil((x_end-x_start)/h)*h+x_start;
            obj.y_end = ceil((y_end-y_start)/h)*h+y_start;
            
            obj.h = h;
            
            obj.n_x = round((obj.x_end-obj.x_start)/h)+1;
            obj.n_y = round((obj.y_end-obj.y_start)/h)+1;  
            
            % reduces round off error
            obj.x_mesh = transpose(round((linspace(obj.x_start, obj.x_end, obj.n_x)-x_start)/h)*h+x_start);
            obj.y_mesh = transpose(round((linspace(obj.y_start, obj.y_end, obj.n_y)-y_start)/h)*h+y_start);
            
            [obj.R_X, obj.R_Y] = meshgrid(obj.x_mesh, obj.y_mesh);
            obj.R_idxs = reshape(1:numel(obj.R_X), size(obj.R_X));

            obj.boundary_X = boundary_X;
            obj.boundary_Y = boundary_Y;

            obj.in_interior = inpolygon_mesh(obj.R_X, obj.R_Y, boundary_X,  boundary_Y);
            
            obj.interior_idxs = obj.R_idxs(obj.in_interior);
            obj.f_R = zeros(obj.n_y, obj.n_x);
        end
        
        function [P, R_patch_idxs] = interpolate_patch(obj, patch, M)
            % constructs vector of idxs of points that are in both patch
            % and cartesian mesh
            [bound_X, bound_Y] = patch.boundary_mesh_xy(false);
            in_patch = inpolygon_mesh(obj.R_X, obj.R_Y, bound_X, bound_Y) & ~obj.in_interior;
            R_patch_idxs = obj.R_idxs(in_patch);
            in_patch = double(in_patch);

            % computing initial "proximity map" with floor and ceil
            % operator
            [XI, ETA] = patch.xi_eta_mesh();
            [patch_X, patch_Y] = patch.convert_to_XY(XI, ETA);

            floor_X_j = floor((patch_X-obj.x_start)/obj.h);
            ceil_X_j = ceil((patch_X-obj.x_start)/obj.h);
            floor_Y_j = floor((patch_Y-obj.y_start)/obj.h);
            ceil_Y_j = ceil((patch_Y-obj.y_start)/obj.h);
            
            n_in_patch = sum(in_patch, 'all');            
            P_xi = zeros(n_in_patch, 1);
            P_eta = zeros(n_in_patch, 1);
            
            P = containers.Map('KeyType', 'int32', 'ValueType', 'any');
            for i = 1:length(R_patch_idxs)
                P_xi(i) = nan; P_eta(i) = nan;
                in_patch(R_patch_idxs(i)) = i;
                P(R_patch_idxs(i)) = nan;
            end
            

            disp("start first pass");
            tic;
            for i = 1:size(patch_X, 1)
                for j = 1:size(patch_X, 2)
                    neighbors = [floor_X_j(i, j), floor_X_j(i, j), ceil_X_j(i, j), ceil_X_j(i, j); floor_Y_j(i, j), ceil_Y_j(i, j), floor_Y_j(i, j), ceil_Y_j(i, j)];
                    
                    for neighbor_i = 1:size(neighbors, 2)
                        neighbor = neighbors(:, neighbor_i) + 1;
                        if any(neighbor > [obj.n_x; obj.n_y]) || any(neighbor < [1; 1])
                            continue;
                        end
                        patch_idx = sub2ind([obj.n_y, obj.n_x], neighbor(2), neighbor(1));
                        
                        if (in_patch(patch_idx) ~= 0) && isnan(P_xi(in_patch(patch_idx)))
                            [xi, eta, converged] = patch.inverse_M_p((neighbor(1)-1)*obj.h + obj.x_start, (neighbor(2)-1)*obj.h + obj.y_start, [XI(i, j); ETA(i, j)]);
                            if converged
                                P_xi(in_patch(patch_idx)) = xi;
                                P_eta(in_patch(patch_idx)) = eta;
                            else
                                warning("Nonconvergence in interpolation")
                            end
                        end
                    end
                end
            end
            toc
            disp("construct nan map")
            
            % second pass for points that aren't touched, could
            % theoretically modify so that we continuously do this until
            % all points are touched
            
            nan_set = containers.Map('KeyType', 'int64', 'ValueType', 'logical');
            % first iterate through and construct nan set            
            for i=1:length(P_xi)
                if isnan(P_xi(i))
                    nan_set(i) = true;
                end
            end
                        
            % pass through nan set until empty
            while nan_set.Count > 0
                for key = keys(nan_set)
                    [i, j] = ind2sub([obj.n_y, obj.n_x], R_patch_idxs(key{1}));
                    
                    neighbor_shifts = [1 -1 0 0 1 1 -1 -1; 0 0 -1 1 1 -1 1 -1];
                    is_touched = false;
                    for neighbor_shift_i = 1:size(neighbor_shifts, 2)
                        neighbor_shift = neighbor_shifts(:, neighbor_shift_i);
                        neighbor = sub2ind([obj.n_y, obj.n_x], i+ neighbor_shift(1), j + neighbor_shift(2));
                        
                        if (in_patch(neighbor) ~= 0) && ~isnan(P_xi(in_patch(neighbor)))
                            [xi, eta, converged] = patch.inverse_M_p((j-1)*obj.h+obj.x_start, (i-1)*obj.h+obj.y_start, [P_xi(in_patch(neighbor)); P_eta(in_patch(neighbor))]);
                            if converged
                                P_xi(key{1}) = xi; P_eta(key{1}) = eta;
                                is_touched = true;
                                break;
                            else
                                warning("Nonconvergence in interpolation")
                            end
                        end
                    end
                    if is_touched
                        remove(nan_set, key{1});
                    end
                end
            end
                        
            f_R_patch = zeros(size(R_patch_idxs));
            for i = 1:length(R_patch_idxs)
                [interior_val, in_range] = patch.locally_compute(P_xi(i), P_eta(i), M);
                if in_range
                     f_R_patch(i) = interior_val;                    
                end
            end
            obj.f_R(R_patch_idxs) = obj.f_R(R_patch_idxs) + f_R_patch;
        end
        
        function fill_interior(obj, f)
        % fills interior of function
            obj.f_R(obj.interior_idxs) = f(obj.R_X(obj.interior_idxs), obj.R_Y(obj.interior_idxs));        
        end
        
        function compute_fc_coeffs(obj)
            obj.fc_coeffs = fftshift(fft2(obj.f_R) / numel(obj.f_R));
        end
        
        function [R_X_err, R_Y_err, f_interpolation, interior_idx] = ifft_interpolation(obj, h_new)
            n_x_err = round((obj.x_end+obj.h-h_new-obj.x_start)/h_new)+1;
            n_y_err = round((obj.y_end+obj.h-h_new-obj.y_start)/h_new)+1;
            
            x_err_mesh = transpose(round((linspace(obj.x_start, obj.x_end+obj.h-h_new, n_x_err)-obj.x_start)/h_new)*h_new+obj.x_start);
            y_err_mesh = transpose(round((linspace(obj.y_start, obj.y_end+obj.h-h_new, n_y_err)-obj.y_start)/h_new)*h_new+obj.y_start);
            
            [R_X_err, R_Y_err] = meshgrid(x_err_mesh, y_err_mesh);

            n_x_diff = n_x_err - obj.n_x;
            n_y_diff = n_y_err - obj.n_y;
            
            padded_fc_coeffs = [zeros(ceil(n_y_diff/2), n_x_err); zeros(size(obj.fc_coeffs, 1), ceil(n_x_diff/2)) obj.fc_coeffs zeros(size(obj.fc_coeffs, 1), floor(n_x_diff/2)); zeros(floor(n_y_diff/2), n_x_err)];
            f_interpolation = real((n_x_err*n_y_err)*ifft2(ifftshift(padded_fc_coeffs)));
            
            idxs = reshape(1:numel(R_X_err), size(R_X_err));
            interior_idx = idxs(inpolygon_mesh(R_X_err, R_Y_err, obj.boundary_X, obj.boundary_Y));
        end
    end
end

