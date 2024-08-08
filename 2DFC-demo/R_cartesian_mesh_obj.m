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
        function obj = R_cartesian_mesh_obj(x_start,x_end, y_start, y_end, h, boundary_X, boundary_Y, ceil_end, in_interior)
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here
            obj.x_start = x_start;
            obj.y_start = y_start;
            if ceil_end
                obj.x_end = ceil((x_end-x_start)/h)*h+x_start;
                obj.y_end = ceil((y_end-y_start)/h)*h+y_start;                
            else
                obj.x_end = round((x_end-x_start)/h)*h+x_start;
                obj.y_end = round((y_end-y_start)/h)*h+y_start;       
            end
            
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
            
            if ~isnan(in_interior)
                if size(in_interior, 1) < obj.n_y || size(in_interior, 2) < obj.n_x
                    in_interior = padarray(in_interior, [obj.n_y - size(in_interior, 1), obj.n_x - size(in_interior, 2)], 'post');
                end
                obj.in_interior = in_interior;
            else
                obj.in_interior = inpolygon_mesh(obj.R_X, obj.R_Y, boundary_X,  boundary_Y);
            end
            
            obj.interior_idxs = obj.R_idxs(obj.in_interior);
            obj.f_R = zeros(obj.n_y, obj.n_x);
        end
        
        function R_new = extend_R(obj, x_start, x_end, y_start, y_end, match_ends)
            if match_ends
                x_start_j = floor((x_start - obj.x_start)/obj.h);
                x_end_j  = ceil((x_end-obj.x_start)/obj.h);

                y_start_j = floor((y_start - obj.y_start)/obj.h);
                y_end_j = ceil((y_end-obj.y_start)/obj.h);
            else
                x_start_j = round((x_start - obj.x_start)/obj.h);
                x_end_j  = round((x_end-obj.x_start)/obj.h);

                y_start_j = round((y_start - obj.y_start)/obj.h);
                y_end_j = round((y_end-obj.y_start)/obj.h);
            end

            matched_in_interior = logical([zeros(-y_start_j, x_end_j-x_start_j+1); zeros(obj.n_y, -x_start_j), obj.in_interior, zeros(obj.n_y, x_end_j-obj.n_x+1); zeros(y_end_j-obj.n_y+1, x_end_j-x_start_j+1)]);

            matched_x_start =x_start_j*obj.h+obj.x_start;
            matched_x_end = x_end_j*obj.h+obj.x_start;
            
            matched_y_start = y_start_j*obj.h+obj.y_start;
            matched_y_end = y_end_j*obj.h+obj.y_start;
                        
            R_new = R_cartesian_mesh_obj(matched_x_start, matched_x_end, matched_y_start, matched_y_end, obj.h, obj.boundary_X, obj.boundary_Y, false, matched_in_interior);
            R_new.f_R(-y_start_j+1:-y_start_j+obj.n_y, -x_start_j+1:-x_start_j+obj.n_x) = obj.f_R;
        end
        
        function [P, R_patch_idxs] = interpolate_patch(obj, patch, M)
            % constructs vector of idxs of points that are in both patch
            % and cartesian mesh
            [bound_X, bound_Y] = patch.boundary_mesh_xy(true);
            in_patch = inpolygon_mesh(obj.R_X, obj.R_Y, bound_X, bound_Y) & ~obj.in_interior;
            R_patch_idxs = obj.R_idxs(in_patch);

            % computing initial "proximity map" with floor and ceil
            % operator
            [XI, ETA] = patch.xi_eta_mesh();
            [patch_X, patch_Y] = patch.convert_to_XY(XI, ETA);

            floor_X_j = floor((patch_X-obj.x_start)/obj.h);
            ceil_X_j = ceil((patch_X-obj.x_start)/obj.h);
            floor_Y_j = floor((patch_Y-obj.y_start)/obj.h);
            ceil_Y_j = ceil((patch_Y-obj.y_start)/obj.h);
            
            P = containers.Map('KeyType', 'int32', 'ValueType', 'any');
            for i = 1:length(R_patch_idxs)
                P(R_patch_idxs(i)) = nan;
            end

            disp("start first pass");

            for i = 1:size(patch_X, 1)
                for j = 1:size(patch_X, 2)
                    neighbors = [floor_X_j(i, j), floor_X_j(i, j), ceil_X_j(i, j), ceil_X_j(i, j); floor_Y_j(i, j), ceil_Y_j(i, j), floor_Y_j(i, j), ceil_Y_j(i, j)];
                    for neighbor_i = 1:size(neighbors, 2)
                        neighbor = neighbors(:, neighbor_i) + 1;
                        if any(neighbor > [obj.n_y; obj.n_x]) || any(neighbor < [1; 1])
                            continue;
                        end
                        patch_idx = sub2ind([obj.n_y, obj.n_x], neighbor(2), neighbor(1));
                        
                        if isKey(P, patch_idx)
                            if isnan(P(patch_idx))
                                [xi, eta, converged] = patch.inverse_M_p((neighbor(1)-1)*obj.h + obj.x_start, (neighbor(2)-1)*obj.h + obj.y_start, [XI(i, j); ETA(i, j)]);
                                if converged
                                    P(patch_idx) = [xi; eta];
                                else
                                    warning("Nonconvergence in interpolation")
                                end
                            end
                        end
                    end
                end
            end
            
            disp("construct nan map")
            
            % second pass for points that aren't touched, could
            % theoretically modify so that we continuously do this until
            % all points are touched
            
            nan_set = containers.Map('KeyType', 'int32', 'ValueType', 'logical');
            % first iterate through and construct nan set            
            for key = keys(P)
                if isnan(P(key{1}))
                    nan_set(key{1}) = true;
                end
            end
                        
            % pass through nan set until empty
            while nan_set.Count > 0
                for key = keys(nan_set)
                    [i, j] = ind2sub([obj.n_y, obj.n_x], key{1});
                    
                    neighbor_shifts = [1 -1 0 0 1 1 -1 -1; 0 0 -1 1 1 -1 1 -1];
                    is_touched = false;
                    for neighbor_shift_i = 1:size(neighbor_shifts, 2)
                        neighbor_shift = neighbor_shifts(:, neighbor_shift_i);
                        neighbor = sub2ind([obj.n_y, obj.n_x], i+ neighbor_shift(1), j + neighbor_shift(2));

                        if isKey(P, neighbor) && ~any(isnan(P(neighbor)))
                            [xi, eta, converged] = patch.inverse_M_p((j-1)*obj.h+obj.x_start, (i-1)*obj.h+obj.y_start, P(neighbor));
                            if converged
                                P(key{1}) = [xi; eta];
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
                xi_eta_point = P(R_patch_idxs(i));
                [interior_val, in_range] = patch.locally_compute(xi_eta_point(1), xi_eta_point(2), M);
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
        
        function [R_X_err, R_Y_err, f_interpolation, interior_idx] = ifft_interpolation(obj, h_new, compute_interior_idxs)
            n_x_err = round((obj.x_end+obj.h-h_new-obj.x_start)/h_new)+1;
            n_y_err = round((obj.y_end+obj.h-h_new-obj.y_start)/h_new)+1;
            
            x_err_mesh = transpose(round((linspace(obj.x_start, obj.x_end+obj.h-h_new, n_x_err)-obj.x_start)/h_new)*h_new+obj.x_start);
            y_err_mesh = transpose(round((linspace(obj.y_start, obj.y_end+obj.h-h_new, n_y_err)-obj.y_start)/h_new)*h_new+obj.y_start);
            
            [R_X_err, R_Y_err] = meshgrid(x_err_mesh, y_err_mesh);

            n_x_diff = n_x_err - obj.n_x;
            n_y_diff = n_y_err - obj.n_y;
            
            padded_fc_coeffs = [zeros(ceil(n_y_diff/2), n_x_err); zeros(size(obj.fc_coeffs, 1), ceil(n_x_diff/2)) obj.fc_coeffs zeros(size(obj.fc_coeffs, 1), floor(n_x_diff/2)); zeros(floor(n_y_diff/2), n_x_err)];
            f_interpolation = real((n_x_err*n_y_err)*ifft2(ifftshift(padded_fc_coeffs)));
            
            if compute_interior_idxs
                idxs = reshape(1:numel(R_X_err), size(R_X_err));
                interior_idx = idxs(inpolygon_mesh(R_X_err, R_Y_err, obj.boundary_X, obj.boundary_Y));
            else
                interior_idx = nan;
            end
        end
    end
end

