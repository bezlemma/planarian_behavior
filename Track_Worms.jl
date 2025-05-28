using TiffImages
using GLMakie
using LinearAlgebra # For norm (distance calculation)
using Statistics    # For mean
using Colors        # For distinguishable_colors, RGB, Gray
using GeometryBasics # For Polygon, Rect

# --- Configuration ---
const THRESHOLD_VALUE = 0.5f0      # Grayscale threshold (0.0 to 1.0)
const MIN_WORM_AREA_PX = 10         # Minimum pixel area to be considered a worm
const MAX_LINKING_DISTANCE_PX = 100.0f0 # Max pixels a worm can move between frames to be linked
const MAX_AVERAGING_DISTANCE_PX = 50.0f0 # Max distance for centroids to be averaged in a frame for single worm assumption

const PIXELS_PER_CM = 58.6f0 # pixels/cm
const MSEC_PER_FRAME = 50.0f0 * 2.0f0 # msec per frame

# --- Helper Functions ---
distance_2d(p1, p2) = norm(p1 - p2)

const CM_PER_PIXEL = 1.0f0 / PIXELS_PER_CM
const AREA_CM2_PER_PIXEL2 = CM_PER_PIXEL^2
const SECONDS_PER_FRAME = MSEC_PER_FRAME / 1000.0f0

struct WormFeatures
    centroid_px::Point2f
    area_px2::Int
    perimeter_px::Float32
    circularity::Float32
    pixels_abs::Vector{Tuple{Int,Int}} # (row, col)
    boundary_pixels_abs::Vector{Point2f} # (col_px, row_px)
end

function find_objects_features(binary_mask::BitMatrix, min_area_px::Int)
    rows, cols = size(binary_mask)
    visited = falses(rows, cols)
    object_features_list = WormFeatures[]

    NEIGHBORS_8 = ((-1, -1), (-1, 0), (-1, 1), (0, -1), (0, 1), (1, -1), (1, 0), (1, 1))
    NEIGHBORS_4 = ((-1,0), (1,0), (0,-1), (0,1))

    for r_start in 1:rows
        for c_start in 1:cols
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
                if obj_area_px2 >= min_area_px
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
                                obj_perimeter_px += 1.0f0
                                is_on_boundary = true
                            end
                        end
                        if is_on_boundary
                             push!(boundary_pixels_for_plot_px, Point2f(pix_c, pix_r))
                        end
                    end
                    obj_circularity = (obj_perimeter_px == 0) ? 0.0f0 : Float32((4 * π * obj_area_px2) / (obj_perimeter_px^2))
                    push!(object_features_list, WormFeatures(Point2f(centroid_x_px, centroid_y_px), obj_area_px2, obj_perimeter_px, obj_circularity, current_object_pixels_tuples, unique(boundary_pixels_for_plot_px))) # Ensure unique boundary pixels
                end
            end
        end
    end
    return object_features_list
end

function find_worm_features_in_frame(img, min_area_px, threshold)
    binary_mask = img .> Float32(threshold)
    return find_objects_features(binary_mask, min_area_px)
end

function average_close_features(features, max_dist)
    if isempty(features)
        return WormFeatures[]
    elseif length(features) == 1
        return features
    end

    avg_centroid_x = mean(f.centroid_px[1] for f in features)
    avg_centroid_y = mean(f.centroid_px[2] for f in features)
    
    all_close = true
    mean_c = Point2f(avg_centroid_x, avg_centroid_y)
    for f in features
        if distance_2d(f.centroid_px, mean_c) > max_dist * 2 
            all_close = false
            break
        end
    end

    if all_close && !isempty(features) 
        avg_area = round(Int, mean(f.area_px2 for f in features))
        avg_perimeter = mean(f.perimeter_px for f in features)
        avg_circularity = mean(f.circularity for f in features)
        
        combined_pixels = vcat((f.pixels_abs for f in features)...)
        combined_boundary_pixels = vcat((f.boundary_pixels_abs for f in features)...) 

        return [WormFeatures(mean_c, avg_area, avg_perimeter, avg_circularity, combined_pixels, unique(combined_boundary_pixels))]
    else 
        return [features[argmax(f.area_px2 for f in features)]]
    end
end


function _core_worm_tracking(frames_iterable, num_frames, img_rows::Int, img_cols::Int; single_worm_assumption::Bool=true)
    last_known_feature = Ref{Union{Nothing, WormFeatures}}(nothing)
    last_known_pos3d_px = Ref{Union{Nothing, Point3f}}(nothing)

    positions_px = Point3f[]
    areas_px2 = Int[]
    perimeters_px = Float32[]
    circularities = Float32[]
    outlines_px = Vector{Point2f}[] 

    for frame_idx in 1:num_frames
        current_frame_img = frames_iterable[frame_idx]
        current_features_list = find_worm_features_in_frame(current_frame_img, MIN_WORM_AREA_PX, THRESHOLD_VALUE)
        
        representative_feature = nothing

        if !isempty(current_features_list)
            if single_worm_assumption
                averaged_features = average_close_features(current_features_list, MAX_AVERAGING_DISTANCE_PX)
                if !isempty(averaged_features)
                    representative_feature = averaged_features[1] 
                end
            else
                representative_feature = current_features_list[argmax(f.area_px2 for f in current_features_list)]
            end
        end

        if !isnothing(representative_feature)
            current_pos_2d_px = representative_feature.centroid_px
            current_pos_3d_px = Point3f(current_pos_2d_px[1], current_pos_2d_px[2], frame_idx)

            if !isnothing(last_known_pos3d_px[])
                dist = distance_2d(Point2f(last_known_pos3d_px[][1:2]), current_pos_2d_px)
                if dist < MAX_LINKING_DISTANCE_PX
                    push!(positions_px, current_pos_3d_px)
                    push!(areas_px2, representative_feature.area_px2)
                    push!(perimeters_px, representative_feature.perimeter_px)
                    push!(circularities, representative_feature.circularity)
                    push!(outlines_px, representative_feature.boundary_pixels_abs)
                    last_known_feature[] = representative_feature
                    last_known_pos3d_px[] = current_pos_3d_px
                else
                    push!(positions_px, current_pos_3d_px) 
                    push!(areas_px2, representative_feature.area_px2)
                    push!(perimeters_px, representative_feature.perimeter_px)
                    push!(circularities, representative_feature.circularity)
                    push!(outlines_px, representative_feature.boundary_pixels_abs)
                    last_known_feature[] = representative_feature
                    last_known_pos3d_px[] = current_pos_3d_px
                end
            else
                push!(positions_px, current_pos_3d_px)
                push!(areas_px2, representative_feature.area_px2)
                push!(perimeters_px, representative_feature.perimeter_px)
                push!(circularities, representative_feature.circularity)
                push!(outlines_px, representative_feature.boundary_pixels_abs)
                last_known_feature[] = representative_feature
                last_known_pos3d_px[] = current_pos_3d_px
            end
        end
        
        if frame_idx % 50 == 0 || frame_idx == num_frames
             println("Core tracking: Processed frame $frame_idx / $num_frames. Found features: $(length(current_features_list))")
        end
    end
    
    return positions_px, areas_px2, perimeters_px, circularities, outlines_px, frames_iterable, img_rows, img_cols
end

# --- Plotting Functions (Separate Windows) ---
function plot_trajectory_3d_single(positions_px_s, img_dims_px, track_color, video_name)
    img_rows, img_cols = img_dims_px
    fig = Figure(size=(1000, 750))
    aspect_xy = (img_rows * CM_PER_PIXEL) / (img_cols * CM_PER_PIXEL) # Simplifies to img_rows/img_cols
    aspect_xy = (img_cols > 0) ? aspect_xy : 1.0 # Avoid division by zero if img_cols is 0

    max_time_val = if !isempty(positions_px_s) 
        maximum(p[3] for p in positions_px_s; init=1.0f0) * SECONDS_PER_FRAME
    else 
        1.0f0 
    end
    max_time_val = (max_time_val == 0.0f0) ? 1.0f0 : max_time_val # ensure not zero
    
    aspect_xz = if (img_cols > 0 && CM_PER_PIXEL > 0)
         max_time_val / (img_cols * CM_PER_PIXEL) 
    else 
        1.0f0 
    end
    aspect_xz = (aspect_xz == 0.0f0) ? 1.0f0 : aspect_xz # ensure not zero

    ax = Axis3(fig[1,1],
               xlabel="X (cm)", ylabel="Y (cm)", zlabel="Time (s)",
               aspect=(1.0, aspect_xy, aspect_xz * 0.3),
               elevation = pi/6)


        positions_cm_s = [Point3f(p[1]*CM_PER_PIXEL, (img_rows - p[2])*CM_PER_PIXEL, p[3]*SECONDS_PER_FRAME) for p in positions_px_s]
        lines!(ax, positions_cm_s, color=track_color, linewidth=3)

    display(GLMakie.Screen(), fig)
    return fig
end

function plot_size_vs_time_single(times_s, areas_px2, track_color, video_name)
    fig = Figure(size=(800, 600))
    ax = Axis(fig[1,1], xlabel="Time (s)", ylabel="Area (cm²)", title="Worm Size Over Time: $video_name")
    if !isempty(times_s)
        areas_cm2 = areas_px2 .* AREA_CM2_PER_PIXEL2
        lines!(ax, times_s, areas_cm2, color=track_color, linewidth=2)
    end
    display(GLMakie.Screen(), fig)
    return fig
end

function plot_circularity_vs_time_single(times_s, circularities, track_color, video_name)
    fig = Figure(size=(800, 600))
    max_circ = isempty(circularities) ? 1.2f0 : maximum(v->isnan(v) ? 0.0f0 : v, circularities; init=1.2f0) * 1.1
    max_circ = max_circ < 0.1f0 ? 1.2f0 : max_circ # Ensure a reasonable upper limit if all values are tiny or NaN
    ax = Axis(fig[1,1], xlabel="Time (s)", ylabel="Circularity (4πA/P²)", title="Worm Circularity Over Time: $video_name", limits=(nothing, nothing, 0, max_circ)) 
    if !isempty(times_s)
        lines!(ax, times_s, circularities, color=track_color, linewidth=2)
    end
    display(GLMakie.Screen(), fig)
    return fig
end


function track_worms_and_features_from_tiff_and_plot_separately(filepath::String; return_raw_positions=false, video_name_suffix="", track_plot_color=RGB(0.1,0.5,0.8))
    println("Loading TIFF stack from: $filepath")
    img_stack_data = TiffImages.load(filepath)

    frames_iterable_local = AbstractMatrix{<:Gray}[]
    num_frames_local = 0
    img_rows_local, img_cols_local = 0, 0

    if ndims(img_stack_data) == 3 && size(img_stack_data, 3) > 1
        num_frames_local = size(img_stack_data, 3)
        img_rows_local, img_cols_local = size(img_stack_data, 1), size(img_stack_data, 2)
        frames_iterable_local = [Gray.(img_stack_data[:, :, i]) for i in 1:num_frames_local]
    elseif isa(img_stack_data, Vector) && !isempty(img_stack_data)
        matrix_elements = filter(x -> isa(x, AbstractMatrix), img_stack_data)
        num_frames_local = length(matrix_elements)
        img_rows_local, img_cols_local = size(matrix_elements[1],1), size(matrix_elements[1],2)
        frames_iterable_local = [Gray.(m) for m in matrix_elements]
    end

    println("Video '$filepath': $num_frames_local frames, $img_rows_local x $img_cols_local px.")

    positions_px, areas_px2, _, circularities, outlines_px, frames_for_outline, r_img, c_img =
        _core_worm_tracking(frames_iterable_local, num_frames_local, img_rows_local, img_cols_local, single_worm_assumption=true)

    if return_raw_positions
        return positions_px, img_rows_local, img_cols_local # MODIFIED: Return image dimensions
    end
    
    video_short_name = basename(filepath) * video_name_suffix


    times_s = [p[3] * SECONDS_PER_FRAME for p in positions_px]

    plot_trajectory_3d_single(positions_px, (r_img, c_img), track_plot_color, video_short_name)
    plot_size_vs_time_single(times_s, areas_px2, track_plot_color, video_short_name)
    plot_circularity_vs_time_single(times_s, circularities, track_plot_color, video_short_name)
        
end


# --- Functions for Video Comparison ---
function plot_path_difference(filepath1::String, filepath2::String)
    println("Calculating path difference between '$filepath1' and '$filepath2'")
    
    raw_data1 = track_worms_and_features_from_tiff_and_plot_separately(filepath1, return_raw_positions=true, video_name_suffix=" (Vid1)")
    raw_data2 = track_worms_and_features_from_tiff_and_plot_separately(filepath2, return_raw_positions=true, video_name_suffix=" (Vid2)")
    
    positions_px_vid1 = raw_data1[1]
    img_rows1, img_cols1 = raw_data1[2], raw_data1[3]
    
    positions_px_vid2 = raw_data2[1]

    println("Video 1 track length (frames): $(length(positions_px_vid1))")
    println("Video 2 track length (frames): $(length(positions_px_vid2))")

    positions_vid1_3d_cm_s = [Point3f(p[1]*CM_PER_PIXEL, (img_rows1 - p[2])*CM_PER_PIXEL, p[3]*SECONDS_PER_FRAME) for p in positions_px_vid1]
    positions_vid2_3d_cm_s = [Point3f(p[1]*CM_PER_PIXEL, (img_rows1 - p[2])*CM_PER_PIXEL, p[3]*SECONDS_PER_FRAME) for p in positions_px_vid2] # Use img_rows1 for y-flip consistency

    map1_time_to_pos_cm = Dict{Float32, Point2f}()
    for p_px in positions_px_vid1 # p_px is (x_px, y_px, frame_idx)
        time_s = p_px[3] * SECONDS_PER_FRAME
        pos_cm = Point2f(p_px[1] * CM_PER_PIXEL, (img_rows1 - p_px[2]) * CM_PER_PIXEL) # Apply y-flip here too for 2D maps
        map1_time_to_pos_cm[time_s] = pos_cm
    end

    map2_time_to_pos_cm = Dict{Float32, Point2f}()
    for p_px in positions_px_vid2
        time_s = p_px[3] * SECONDS_PER_FRAME
        pos_cm = Point2f(p_px[1] * CM_PER_PIXEL, (img_rows1 - p_px[2]) * CM_PER_PIXEL) # Apply y-flip
        map2_time_to_pos_cm[time_s] = pos_cm
    end
    
    common_times_s = sort(collect(intersect(keys(map1_time_to_pos_cm), keys(map2_time_to_pos_cm))))

    println("Found $(length(common_times_s)) common time points for difference calculation.")

    differences_cm_magnitude = Float32[]
    difference_track_3d_points = Point3f[]

    for t_s in common_times_s
        pos1_2d_cm = map1_time_to_pos_cm[t_s]
        pos2_2d_cm = map2_time_to_pos_cm[t_s]
        
        diff_vec_2d_cm = pos1_2d_cm - pos2_2d_cm
        push!(differences_cm_magnitude, norm(diff_vec_2d_cm))
        push!(difference_track_3d_points, Point3f(diff_vec_2d_cm[1], diff_vec_2d_cm[2], t_s))
    end

    # --- Plotting ---
    fig_diff = Figure(size=(1200, 1000))
    
    max_time_s_orig_1 = isempty(positions_vid1_3d_cm_s) ? 1.0f0 : maximum(p[3] for p in positions_vid1_3d_cm_s; init=1.0f0)
    max_time_s_orig_2 = isempty(positions_vid2_3d_cm_s) ? 1.0f0 : maximum(p[3] for p in positions_vid2_3d_cm_s; init=1.0f0)
    max_time_s_orig = max(max_time_s_orig_1, max_time_s_orig_2, 1.0f0)
    max_time_s_orig = (max_time_s_orig == 0.0f0) ? 1.0f0 : max_time_s_orig
    
    aspect_xz_orig = (img_cols1 > 0 && CM_PER_PIXEL > 0) ? max_time_s_orig / (img_cols1 * CM_PER_PIXEL) : 1.0
    aspect_xz_orig = (aspect_xz_orig == 0.0f0) ? 1.0f0 : aspect_xz_orig


    ax_tracks_3d = Axis3(fig_diff[1,1], title="Original Tracks 3D",
                         xlabel="X (cm)", ylabel="Y (cm)", zlabel="Time (s)")
    lines!(ax_tracks_3d, positions_vid1_3d_cm_s, color=:blue, linewidth=2)
    lines!(ax_tracks_3d, positions_vid2_3d_cm_s, color=:red, linewidth=2)

    max_abs_dx = isempty(difference_track_3d_points) ? 1.0f0 : maximum(abs(p[1]) for p in difference_track_3d_points; init=1.0f0)
    max_abs_dy = isempty(difference_track_3d_points) ? 1.0f0 : maximum(abs(p[2]) for p in difference_track_3d_points; init=1.0f0)
    max_time_diff = isempty(difference_track_3d_points) ? 1.0f0 : maximum(p[3] for p in difference_track_3d_points; init=1.0f0)

    max_abs_dx = (max_abs_dx == 0.0f0) ? 1.0f0 : max_abs_dx # Avoid division by zero for aspect
    max_abs_dy = (max_abs_dy == 0.0f0) ? 1.0f0 : max_abs_dy
    max_time_diff = (max_time_diff == 0.0f0) ? 1.0f0 : max_time_diff

    ax_diff_track_3d = Axis3(fig_diff[1,2], xlabel="dX (cm)", ylabel="dY (cm)", zlabel="Time (s)", title="Trail Diff")
    lines!(ax_diff_track_3d, difference_track_3d_points, color=:green, linewidth=2)
    
    # 3. Plot absolute magnitude of 2D difference
    ax_diff_val = Axis(fig_diff[2,1:2], xlabel="Time (s)", ylabel="Path Difference Magnitude (cm)", title="Absolute Magnitude of Diff")
    lines!(ax_diff_val, common_times_s, differences_cm_magnitude, color=:purple, linewidth=2.5)
    
    display(GLMakie.Screen(), fig_diff)
    return fig_diff
end

TIFF_FILE_PATH_VID1 = "/Users/dl0346/Documents/PlanarianVideos/ColloidTestInverse.tif" 
TIFF_FILE_PATH_VID2 = "/Users/dl0346/Documents/PlanarianVideos/ColloidInverse2.tif" 

 plot_path_difference(TIFF_FILE_PATH_VID1, TIFF_FILE_PATH_VID2)