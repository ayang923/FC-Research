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
        f_R
        
        fc_coeffs
    end
    
    methods
        function obj = R_cartesian_mesh_obj(x_start,x_end, y_start, y_end, h, boundary_X, boundary_Y)
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here
            obj.x_start = x_start;
            obj.x_end = ceil((x_end-x_start)/h)*h+x_start;
            obj.y_start = y_start;
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
            obj.in_interior = inpolygon(obj.R_X, obj.R_Y, boundary_X,  boundary_Y);
            obj.f_R = zeros(obj.n_y, obj.n_x);
        end
        
        function interpolate_patch(obj, patch, d)
            hash_sig = 14;
            % constructs vector of idxs of points that are in both patch
            % and cartesian mesh
            [bound_X, bound_Y] = patch.boundary_mesh_xy();
            in_patch = inpolygon(obj.R_X, obj.R_Y, bound_X, bound_Y) & ~obj.in_interior;
            R_patch_idxs = obj.R_idxs(in_patch);
            
            R_patch_X = obj.R_X(R_patch_idxs);
            R_patch_Y = obj.R_Y(R_patch_idxs);            
            
            % computing initial "proximity map" with floor and ceil
            % operator
            [XI, ETA] = patch.xi_eta_mesh();
            [patch_X, patch_Y] = patch.convert_to_XY(XI, ETA);

            floor_X = floor((patch_X-obj.x_start)/obj.h)*obj.h + obj.x_start;
            ceil_X = ceil((patch_X-obj.x_start)/obj.h)*obj.h + obj.x_start;
            floor_Y = floor((patch_Y-obj.y_start)/obj.h)*obj.h + obj.y_start;
            ceil_Y = ceil((patch_Y-obj.y_start)/obj.h)*obj.h + obj.y_start;
                        
            P = containers.Map('KeyType', 'char', 'ValueType', 'any');
            for i = 1:length(R_patch_X)
                P(mat2str([R_patch_X(i); R_patch_Y(i)], hash_sig)) = nan;
            end

            for i = 1:size(patch_X, 1)
                for j = 1:size(patch_X, 2)
                    neighbors = [floor_X(i, j), floor_X(i, j), ceil_X(i, j), ceil_X(i, j); floor_Y(i, j), ceil_Y(i, j), floor_Y(i, j), ceil_Y(i, j)];
                    for neighbor_i = 1:size(neighbors, 2)
                        neighbor = neighbors(:, neighbor_i);
                        str_neighbor = mat2str(neighbor, hash_sig);
                        
                        if isKey(P, str_neighbor)
                            if isnan(P(str_neighbor))
                                [xi, eta, converged] = patch.inverse_M_p(neighbor(1), neighbor(2), [XI(i, j); ETA(i, j)]);
                                if converged
                                    P(str_neighbor) = [xi; eta];
                                else
                                    warning("Nonconvergence in interpolation")
                                end
                            end
                        end
                    end
                end
            end          
            
            % second pass for points that aren't touched, could
            % theoretically modify so that we continuously do this until
            % all points are touched
            nan_set = containers.Map('KeyType', 'char', 'ValueType', 'logical');
            while true
                for key = keys(P)
                    if isnan(P(key{1}))
                        is_touched = false;
                        
                        pnt = eval(key{1});
                        neighboridxs = [1 -1 0 0 1 1 -1 -1; 0 0 -1 1 1 -1 1 -1];
                        for neighboridx_i = 1:size(neighboridxs, 2)
                            neighboridx = neighboridxs(:, neighboridx_i);
                            neighbor = pnt + neighboridx * obj.h;
                            str_neighbor = mat2str(neighbor, hash_sig);
                            if isKey(P, str_neighbor) && ~any(isnan((P(str_neighbor))))
                                [xi, eta, converged] = patch.inverse_M_p(pnt(1), pnt(2), P(str_neighbor));
                                if converged
                                    P(key{1}) = [xi; eta];

                                    is_touched = true;
                                    if isKey(nan_set, key{1})
                                        remove(nan_set, key{1});
                                    end
                                    break;
                                else
                                    warning("Nonconvergence in interpolation")
                                end
                            end
                        end
                        if ~is_touched
                            nan_set(key{1}) = true;
                        end
                    end
                end
                if isempty(keys(nan_set))
                    break;
                end
            end
            
            f_R_patch = zeros(size(R_patch_X));
            for i = 1:length(R_patch_X)
                xi_eta_point = P(mat2str([R_patch_X(i); R_patch_Y(i)], hash_sig));
                [f_R_patch(i), in_range] = patch.locally_compute(xi_eta_point(1), xi_eta_point(2), d);
                if ~in_range
                    disp('huh')
                end
            end
            obj.f_R(R_patch_idxs) = obj.f_R(R_patch_idxs) + f_R_patch;
        end
        
        function fill_interior(obj, f)
        % fills interior of function
            interior_idxs = obj.R_idxs(obj.in_interior);
            obj.f_R(interior_idxs) = f(obj.R_X(interior_idxs), obj.R_Y(interior_idxs));        
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
            f_interpolation = real((n_x_err*n_y_err)*ifft2(fftshift(padded_fc_coeffs)));
            
            idxs = reshape(1:numel(R_X_err), size(R_X_err));
            interior_idx = idxs(inpolygon(R_X_err, R_Y_err, obj.boundary_X, obj.boundary_Y));
        end
    end
end

