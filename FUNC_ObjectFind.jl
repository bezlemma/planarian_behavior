# Module to detect connected “worm” objects in a binary image and compute their geometric features.
module FUNC_ObjectFind

using LinearAlgebra, Statistics, GeometryBasics 

export find_objects_features, remove_outlier_objects

struct WormFeaturesLocal
    centroid::Point2f
    wormarea::Int
    major_axis::Float32
    minor_axis::Float32
    pixels_abs::Vector{Tuple{Int,Int}}
    boundary_pixels_abs::Vector{Point2f}
end

function find_objects_features(binary_mask, MIN_AREA, MAX_AREA, MAJOR_MINOR_MAX_RATIO)

    rows, cols = size(binary_mask)     # get image dimensions

    visited = falses(rows, cols)     # track visited pixels

    object_features_list = WormFeaturesLocal[]    # output list
    
    # neighbor offsets for flood‐fill (8‐connectivity) and perimeter check (4)
    NEIGHBORS_8 = ((-1, -1), (-1, 0), (-1, 1),
                   ( 0, -1),          ( 0, 1),
                   ( 1, -1), ( 1, 0), ( 1, 1))
    NEIGHBORS_4 = ((-1,0), (1,0), (0,-1), (0,1))

    # scan every pixel
    for r_start in 1:rows, c_start in 1:cols

        # start a flood‐fill if pixel is object and unvisited
        if binary_mask[r_start, c_start] && !visited[r_start, c_start]
            current_object = Tuple{Int,Int}[]
            q = Tuple{Int,Int}[(r_start, c_start)]
            visited[r_start, c_start] = true
            head = 1

            # BFS to collect all connected pixels
            while head <= length(q)
                curr_r, curr_c = q[head]; head += 1
                push!(current_object, (curr_r, curr_c))
                for (dr, dc) in NEIGHBORS_8
                    nr, nc = curr_r + dr, curr_c + dc
                    if (1 <= nr <= rows) && (1 <= nc <= cols) &&
                       binary_mask[nr, nc] && !visited[nr, nc]
                        visited[nr, nc] = true
                        push!(q, (nr, nc))
                    end
                end
            end

            # compute object area
            obj_area = length(current_object)

            # filter by size thresholds
            if obj_area >= MIN_AREA && obj_area <= MAX_AREA
                
                # compute centroid by averaging row/col
                sum_x, sum_y = 0.0, 0.0
                for (pix_r, pix_c) in current_object
                    sum_x += pix_c; sum_y += pix_r
                end
                centroid_x = sum_x / obj_area
                centroid_y = sum_y / obj_area

                # compute covariance matrix for ellipse fitting
                xs = [pix_c for (pix_r, pix_c) in current_object]
                ys = [pix_r for (pix_r, pix_c) in current_object]
                cx, cy = centroid_x, centroid_y
                cov_xx = mean((xs .- cx).^2)
                cov_yy = mean((ys .- cy).^2)
                cov_xy = mean((xs .- cx) .* (ys .- cy))
                M = Symmetric(Float32[cov_xx cov_xy; cov_xy cov_yy])
                ev = eigen(M)
                λ = sort(ev.values; rev=true)
                major_axis = 2 * sqrt(λ[1])
                minor_axis = 2 * sqrt(λ[2])

                # filter by major/minor ratio
                if (major_axis/major_axis) > MAJOR_MINOR_MAX_RATIO
                    continue
                end

                # compute boundary pixels for plotting
                obj_perimeter = 0.0
                boundary_pixels = Point2f[]
                pixel_set = Set(current_object)
                for (pix_r, pix_c) in current_object
                    is_on_boundary = false
                    for (dr, dc) in NEIGHBORS_4
                        nr_edge, nc_edge = pix_r + dr, pix_c + dc
                        if !((1 <= nr_edge <= rows) && (1 <= nc_edge <= cols) && ((nr_edge, nc_edge) in pixel_set))
                            obj_perimeter += 1.0
                            is_on_boundary = true
                        end
                    end
                    if is_on_boundary
                        push!(boundary_pixels, Point2f(pix_c, pix_r))
                    end
                end

                # save the feature struct
                push!(object_features_list,
                      WormFeaturesLocal(Point2f(centroid_x, centroid_y),
                                        obj_area, major_axis, minor_axis,
                                        current_object,
                                        unique(boundary_pixels)))
            end
        end
    end

    # return all detected worm features
    return object_features_list
end

# Allow calling with 3 args, using no axis filtering
function find_objects_features(binary_mask, MIN_AREA, MAX_AREA)
    return find_objects_features(binary_mask, MIN_AREA, MAX_AREA, typemax(Float32))
end
# Allow calling with 1 args, using no filtering
function find_objects_features(binary_mask)
    return find_objects_features(binary_mask, 0.0, typemax(Float32), typemax(Float32))
end

#remove objects whose area is outside [MIN_AREA, MAX_AREA]
function remove_outlier_objects(binary_mask::BitArray{3}, MIN_AREA, MAX_AREA, MAJOR_MINOR_MAX_RATIO)
    rows, cols, slices = size(binary_mask)
    filtered_mask = falses(rows, cols, slices)

    for s in 1:slices
        # Process each 2D slice
        current_slice_mask = binary_mask[:, :, s]
        # find all features within area bounds for the current slice
        valid_feats = find_objects_features(current_slice_mask, MIN_AREA, MAX_AREA, MAJOR_MINOR_MAX_RATIO)
        
        # reconstruct mask from valid object pixels for the current slice
        for feat in valid_feats
            for (r, c) in feat.pixels_abs
                filtered_mask[r, c, s] = true
            end
        end
    end
    return filtered_mask
end


end # module