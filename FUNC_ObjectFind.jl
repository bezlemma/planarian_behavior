# Module to detect connected “worm” objects in a binary image
# and compute their geometric features.
module FUNC_ObjectFind

using LinearAlgebra
using Statistics    
using GeometryBasics 

export find_objects_features

struct WormFeaturesLocal
    centroid::Point2f
    wormarea::Int
    major_axis::Float32
    minor_axis::Float32
    pixels_abs::Vector{Tuple{Int,Int}}
    boundary_pixels_abs::Vector{Point2f}
end

function find_objects_features(binary_mask, min_area_px, MAX_AREA_PX, min_major_axis, max_major_axis, min_minor_axis, max_minor_axis)

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
            current_object_pixels_tuples = Tuple{Int,Int}[]
            q = Tuple{Int,Int}[(r_start, c_start)]
            visited[r_start, c_start] = true
            head = 1

            # BFS to collect all connected pixels
            while head <= length(q)
                curr_r, curr_c = q[head]; head += 1
                push!(current_object_pixels_tuples, (curr_r, curr_c))
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
            obj_area_px2 = length(current_object_pixels_tuples)

            # filter by size thresholds
            if obj_area_px2 >= min_area_px && obj_area_px2 <= MAX_AREA_PX
                
                # compute centroid by averaging row/col
                sum_x, sum_y = 0.0f0, 0.0f0
                for (pix_r, pix_c) in current_object_pixels_tuples
                    sum_x += pix_c; sum_y += pix_r
                end
                centroid_x_px = sum_x / obj_area_px2
                centroid_y_px = sum_y / obj_area_px2

                # compute covariance matrix for ellipse fitting
                xs = [pix_c for (pix_r, pix_c) in current_object_pixels_tuples]
                ys = [pix_r for (pix_r, pix_c) in current_object_pixels_tuples]
                cx, cy = centroid_x_px, centroid_y_px
                cov_xx = mean((xs .- cx).^2)
                cov_yy = mean((ys .- cy).^2)
                cov_xy = mean((xs .- cx) .* (ys .- cy))
                M = Symmetric(Float32[cov_xx cov_xy; cov_xy cov_yy])
                ev = eigen(M)
                λ = sort(ev.values; rev=true)
                major_axis = 2f0 * sqrt(λ[1])
                minor_axis = 2f0 * sqrt(λ[2])

                # filter by major/minor axis thresholds
                if !(major_axis >= min_major_axis && major_axis <= max_major_axis &&
                     minor_axis >= min_minor_axis && minor_axis <= max_minor_axis)
                    continue
                end

                # compute boundary pixels for plotting (optional)
                obj_perimeter_px = 0.0f0
                boundary_pixels_for_plot_px = Point2f[]
                pixel_set = Set(current_object_pixels_tuples)
                for (pix_r, pix_c) in current_object_pixels_tuples
                    is_on_boundary = false
                    for (dr, dc) in NEIGHBORS_4
                        nr_edge, nc_edge = pix_r + dr, pix_c + dc
                        if !((1 <= nr_edge <= rows) && (1 <= nc_edge <= cols) && ((nr_edge, nc_edge) in pixel_set))
                            obj_perimeter_px += 1.0f0
                            is_on_boundary = true
                        end
                    end
                    if is_on_boundary
                        push!(boundary_pixels_for_plot_px, Point2f(pix_c, pix_r))
                    end
                end

                # save the feature struct
                push!(object_features_list,
                      WormFeaturesLocal(Point2f(centroid_x_px, centroid_y_px),
                                        obj_area_px2, major_axis, minor_axis,
                                        current_object_pixels_tuples,
                                        unique(boundary_pixels_for_plot_px)))
            end
        end
    end

    # return all detected worm features
    return object_features_list
end

# Backwards compatibility: allow calling with 3 args, using no axis filtering
function find_objects_features(binary_mask, min_area_px, MAX_AREA_PX)
    return find_objects_features(binary_mask, min_area_px, MAX_AREA_PX,
                                 0.0f0, typemax(Float32),
                                 0.0f0, typemax(Float32))
end
end # module