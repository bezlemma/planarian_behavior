module FUNC_ObjectFind

using LinearAlgebra # For norm (distance calculation)
using Statistics    # For mean
using GeometryBasics # For Polygon, Rect, Point2f, Point3f, Circle\

export find_objects_features

struct WormFeaturesLocal
    centroid::Point2f
    wormarea::Int
    wormperimeter::Float32
    circularity::Float32
    pixels_abs::Vector{Tuple{Int,Int}} # (row, col)
    boundary_pixels_abs::Vector{Point2f}
end

function find_objects_features(binary_mask, min_area_px, MAX_AREA_PX) 
    rows, cols = size(binary_mask)
    visited = falses(rows, cols)
    object_features_list = WormFeaturesLocal[]
    NEIGHBORS_8 = ((-1, -1), (-1, 0), (-1, 1), (0, -1), (0, 1), (1, -1), (1, 0), (1, 1))
    NEIGHBORS_4 = ((-1,0), (1,0), (0,-1), (0,1))

    for r_start in 1:rows, c_start in 1:cols
        if binary_mask[r_start, c_start] && !visited[r_start, c_start]
            current_object_pixels_tuples = Tuple{Int,Int}[]
            q = Tuple{Int,Int}[(r_start, c_start)]
            visited[r_start, c_start] = true
            head = 1
            while head <= length(q)
                curr_r, curr_c = q[head]; head += 1
                push!(current_object_pixels_tuples, (curr_r, curr_c))
                for (dr, dc) in NEIGHBORS_8
                    nr, nc = curr_r + dr, curr_c + dc
                    if (1 <= nr <= rows) && (1 <= nc <= cols) && binary_mask[nr, nc] && !visited[nr, nc]
                        visited[nr, nc] = true
                        push!(q, (nr, nc))
                    end
                end
            end
            obj_area_px2 = length(current_object_pixels_tuples)

            if obj_area_px2 >= min_area_px && obj_area_px2 <= MAX_AREA_PX
                sum_x, sum_y = 0.0f0, 0.0f0
                for (pix_r, pix_c) in current_object_pixels_tuples
                    sum_x += pix_c; sum_y += pix_r
                end
                centroid_x_px = sum_x / obj_area_px2
                centroid_y_px = sum_y / obj_area_px2
                obj_perimeter_px = 0.0f0
                boundary_pixels_for_plot_px = Point2f[]
                pixel_set = Set(current_object_pixels_tuples)
                for (pix_r, pix_c) in current_object_pixels_tuples
                    is_on_boundary = false
                    for (dr, dc) in NEIGHBORS_4
                        nr_edge, nc_edge = pix_r + dr, pix_c + dc
                        if !((1 <= nr_edge <= rows) && (1 <= nc_edge <= cols) && ((nr_edge, nc_edge) in pixel_set))
                            obj_perimeter_px += 1.0f0; is_on_boundary = true
                        end
                    end
                    if is_on_boundary push!(boundary_pixels_for_plot_px, Point2f(pix_c, pix_r)) end
                end
                obj_circularity = (obj_perimeter_px == 0) ? 0.0f0 : Float32((4 * π * obj_area_px2) / (obj_perimeter_px^2))
                push!(object_features_list, WormFeaturesLocal(Point2f(centroid_x_px, centroid_y_px), obj_area_px2, obj_perimeter_px, obj_circularity, current_object_pixels_tuples, unique(boundary_pixels_for_plot_px)))
            end
        end
    end
    return object_features_list
end


end # module