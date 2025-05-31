using TiffImages
using GLMakie
using LinearAlgebra # For norm (distance calculation)
using Statistics    # For mean
using Colors        # For distinguishable_colors, RGB, Gray
using GeometryBasics # For Polygon, Rect, Point2f, Point3f, Circle
using AbstractFFTs # For FFT-based phase correlation
using FFTW       # FFTW is often the backend for AbstractFFTs
using ProgressMeter # For progress bar
using Base.Threads  # For multithreading
using Observables   # For interactive Observables with GLMakie

# --- Constants ---
const THRESHOLD_VALUE = 1.5f0
const MIN_WORM_AREA_PX = 15
const WORM_DISTANCE_LINK = 40.0f0       # Base linking distance
const WORM_DISTANCE_SEARCHADD = 0.1f0  # Expansion per missed frame for search radius
const MIN_ACCEPTABLE_CIRCULARITY = 0.10f0
const MAX_AREA_PX = 500 # NEW: Maximum acceptable worm area to filter out large plate artifacts
const MIN_AREA_PX = MIN_WORM_AREA_PX # Redefine for clarity below if needed, but MIN_WORM_AREA_PX is used.

const CM_PER_PIXEL = 1.0f0 / 58.6f0
const MSEC_PER_FRAME = 50.0f0 * 2.0f0

# --- Structs ---
struct WormFeatures
    centroid_px::Point2f
    area_px2::Int
    perimeter_px::Float32
    circularity::Float32
    pixels_abs::Vector{Tuple{Int,Int}} # (row, col)
    boundary_pixels_abs::Vector{Point2f} # (col_px, row_px)
end

struct WormProcessedData
    filepath::String
    name::String
    positions_px::Vector{Point3f} # (x_px, y_px, frame_idx) - frame_idx is the actual frame number
    areas_px2::Vector{Int}
    perimeters_px::Vector{Float32}
    circularities::Vector{Float32}
    img_rows::Int
    img_cols::Int
    times_s::Vector{Float32}
    positions_cm::Vector{Point3f} # (x_cm, y_cm_flipped, time_s)
end

# --- Helper Functions (Original) ---
distance_2d(p1::Point2f, p2::Point2f) = norm(p1 - p2)

function find_objects_features(binary_mask, min_area_px) # Removed ::BitMatrix for flexibility if needed
    rows, cols = size(binary_mask)
    visited = falses(rows, cols)
    object_features_list = WormFeatures[]
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
            # --- MODIFICATION: Add MAX_AREA_PX check ---
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
                push!(object_features_list, WormFeatures(Point2f(centroid_x_px, centroid_y_px), obj_area_px2, obj_perimeter_px, obj_circularity, current_object_pixels_tuples, unique(boundary_pixels_for_plot_px)))
            end
        end
    end
    return object_features_list
end

# --- Preprocessing Helpers (Unchanged) ---
# ... (create_border_image, fftshift2d, phase_correlate, translate_image remain the same)
function create_border_image(img::AbstractMatrix{T}) where T
    rows, cols = size(img)
    border_only_img = copy(img)
    center_h = floor(Int, rows * 0.9)
    center_w = floor(Int, cols * 0.9)
    if center_h < rows && center_w < cols && center_h > 0 && center_w > 0
        margin_y = rows - center_h; offset_y = margin_y ÷ 2; start_y_center = offset_y + 1
        end_y_center = start_y_center + center_h - 1
        margin_x = cols - center_w; offset_x = margin_x ÷ 2; start_x_center = offset_x + 1
        end_x_center = start_x_center + center_w - 1
        actual_start_y = max(1, start_y_center); actual_end_y = min(rows, end_y_center)
        actual_start_x = max(1, start_x_center); actual_end_x = min(cols, end_x_center)
        if actual_start_y <= actual_end_y && actual_start_x <= actual_end_x
            border_only_img[actual_start_y:actual_end_y, actual_start_x:actual_end_x] .= T(0)
        end
    elseif center_h >= rows || center_w >= cols
        border_only_img .= T(0)
    end
    return border_only_img
end
function fftshift2d(mat) sy, sx = size(mat); return circshift(mat, (sy ÷ 2, sx ÷ 2)) end
function phase_correlate(ref_img_orig, moving_img_orig)
    ref_img_cropped, moving_img_cropped = if size(ref_img_orig) != size(moving_img_orig)
        r_h, r_w = size(ref_img_orig); m_h, m_w = size(moving_img_orig)
        min_h = min(r_h, m_h); min_w = min(r_w, m_w)
        ref_img_orig[1:min_h, 1:min_w], moving_img_orig[1:min_h, 1:min_w]
    else ref_img_orig, moving_img_orig end
    ref_img_border = create_border_image(ref_img_cropped)
    moving_img_border = create_border_image(moving_img_cropped)
    fft_ref = fft(Float32.(ref_img_border)); fft_moving = fft(Float32.(moving_img_border))
    cross_power_spectrum = fft_ref .* conj(fft_moving)
    magnitude = abs.(cross_power_spectrum)
    epsilon = Float32(1e-10)
    normalized_cross_power_spectrum = cross_power_spectrum ./ (magnitude .+ epsilon)
    correlation_matrix = fftshift2d(real.(ifft(normalized_cross_power_spectrum)))
    _max_val, max_idx_linear = findmax(correlation_matrix)
    max_idx_cartesian = CartesianIndices(size(correlation_matrix))[max_idx_linear]
    center_y, center_x = (size(correlation_matrix) .÷ 2) .+ 1
    dy = max_idx_cartesian[1] - center_y; dx = max_idx_cartesian[2] - center_x
    return dy, dx
end
function translate_image(img::AbstractMatrix{T}, dy, dx, fill_value::T) where T
    rows, cols = size(img); translated_img = fill(fill_value, rows, cols)
    src_y_start = max(1, 1 - dy); src_y_end = min(rows, rows - dy)
    src_x_start = max(1, 1 - dx); src_x_end = min(cols, cols - dx)
    dest_y_start = max(1, 1 + dy); dest_y_end = min(rows, rows + dy)
    dest_x_start = max(1, 1 + dx); dest_x_end = min(cols, cols + dx)
    if (src_y_start <= src_y_end) && (src_x_start <= src_x_end) && (dest_y_start <= dest_y_end) && (dest_x_start <= dest_x_end)
        src_height = src_y_end - src_y_start + 1; src_width = src_x_end - src_x_start + 1
        dest_height = dest_y_end - dest_y_start + 1; dest_width = dest_x_end - dest_x_start + 1
        copy_height = min(src_height, dest_height); copy_width = min(src_width, dest_width)
        src_y_end_clipped = src_y_start + copy_height - 1; src_x_end_clipped = src_x_start + copy_width - 1
        dest_y_end_clipped = dest_y_start + copy_height - 1; dest_x_end_clipped = dest_x_start + copy_width - 1
        if copy_height > 0 && copy_width > 0
            translated_img[dest_y_start:dest_y_end_clipped, dest_x_start:dest_x_end_clipped] = img[src_y_start:src_y_end_clipped, src_x_start:src_x_end_clipped]
        end
    end
    return translated_img
end

# --- Circular Interpolation (Unchanged) ---
# ... (interpolate_circular_path remains the same)
function interpolate_circular_path(p_start_px::Point2f, p_end_px::Point2f, circle_center_px::Point2f, circle_radius_px::Float32, num_intermediate_points::Int)::Vector{Point2f}
    if num_intermediate_points <= 0 return Point2f[] end
    v_start = p_start_px - circle_center_px; v_end = p_end_px - circle_center_px
    angle_start_rad = atan(v_start[2], v_start[1]); angle_end_rad = atan(v_end[2], v_end[1])
    delta_angle_rad = angle_end_rad - angle_start_rad
    if delta_angle_rad > π delta_angle_rad -= 2π elseif delta_angle_rad < -π delta_angle_rad += 2π end
    interpolated_points = Vector{Point2f}(undef, num_intermediate_points)
    total_steps_for_interpolation = num_intermediate_points + 1
    for i in 1:num_intermediate_points
        fraction = Float32(i) / Float32(total_steps_for_interpolation)
        current_angle_rad = angle_start_rad + fraction * delta_angle_rad
        pt_x = circle_center_px[1] + circle_radius_px * cos(current_angle_rad)
        pt_y = circle_center_px[2] + circle_radius_px * sin(current_angle_rad)
        interpolated_points[i] = Point2f(pt_x, pt_y)
    end
    return interpolated_points
end

# --- Function to generate the full binary stack (Unchanged) ---
# ... (generate_full_binary_stack remains the same)
function generate_full_binary_stack(float_stack::AbstractArray{Float32, 3}, threshold::Float32)::BitArray{3}
    rows, cols, num_frames = size(float_stack)
    binary_stack = BitArray{3}(undef, rows, cols, num_frames)
    p_binarize = Progress(num_frames, "Binarizing full stack: ")
    for i in 1:num_frames
        binary_stack[:, :, i] = float_stack[:, :, i] .> threshold
        next!(p_binarize)
    end
    finish!(p_binarize)
    return binary_stack
end

# --- Core Tracking Logic (MODIFIED) ---
function _core_worm_tracking_single_pass(
    binary_stack_global::BitArray{3},
    num_total_frames::Int,
    img_rows::Int, img_cols::Int,
    is_forward_pass::Bool,
    img_center_px::Point2f, arena_radius_px::Float32,
    on_circle_tolerance_px::Float32, arena_boundary_tolerance_px::Float32
)::Vector{Union{Nothing,WormFeatures}}

    pass_features_output = Vector{Union{Nothing, WormFeatures}}(nothing, num_total_frames)
    _last_identified_feature_this_pass::Union{Nothing, WormFeatures} = nothing
    _last_pos_for_search_expansion_this_pass::Union{Nothing, Point2f} = nothing
    _consecutive_misses_this_pass::Int = 0
    frame_iterator = is_forward_pass ? (1:num_total_frames) : (num_total_frames:-1:1)

    for current_frame_idx in frame_iterator
        binary_mask_this_frame = view(binary_stack_global, :, :, current_frame_idx)
        # find_objects_features now also filters by MAX_AREA_PX
        all_objects_in_frame = find_objects_features(binary_mask_this_frame, MIN_WORM_AREA_PX)
        
        found_feature_this_frame::Union{Nothing, WormFeatures} = nothing
        candidate_features_this_frame = WormFeatures[]
        search_center::Point2f = Point2f(0,0) # Ensure search_center is defined

        if !isnothing(_last_identified_feature_this_pass)
            search_center = _last_identified_feature_this_pass.centroid_px
            search_radius = WORM_DISTANCE_LINK 
            candidate_features_this_frame = filter(
                f -> distance_2d(f.centroid_px, search_center) <= search_radius &&
                     distance_2d(f.centroid_px, img_center_px) <= arena_radius_px + arena_boundary_tolerance_px,
                all_objects_in_frame)
        elseif !isnothing(_last_pos_for_search_expansion_this_pass)
            search_center = _last_pos_for_search_expansion_this_pass
            search_radius = WORM_DISTANCE_LINK + _consecutive_misses_this_pass * WORM_DISTANCE_SEARCHADD
            potential_candidates = filter(
                f -> distance_2d(f.centroid_px, search_center) <= search_radius &&
                     distance_2d(f.centroid_px, img_center_px) <= arena_radius_px + arena_boundary_tolerance_px,
                all_objects_in_frame)
            # When "lost", prioritize objects near the circular arena edge
            # This might be too restrictive if worms are lost in the center. Consider alternatives if this fails.
            candidate_features_this_frame = filter(
                f -> abs(distance_2d(f.centroid_px, img_center_px) - arena_radius_px) <= on_circle_tolerance_px,
                potential_candidates)
            # If the above yields no candidates (e.g. worm lost in center), fall back to `potential_candidates`?
            # For now, keeping the strict edge search when lost.
        else 
            search_center = img_center_px # search_center used for scoring if multiple candidates
            search_radius = arena_radius_px + arena_boundary_tolerance_px
            candidate_features_this_frame = filter(
                f -> distance_2d(f.centroid_px, img_center_px) <= search_radius,
                all_objects_in_frame)
        end

        if !isempty(candidate_features_this_frame)
            if length(candidate_features_this_frame) == 1
                found_feature_this_frame = candidate_features_this_frame[1]
            else
                best_score = Inf32
                best_candidate_idx = -1

                norm_dist_factor = WORM_DISTANCE_LINK # Default normalization for distance
                if isnothing(_last_identified_feature_this_pass) && isnothing(_last_pos_for_search_expansion_this_pass)
                    norm_dist_factor = search_radius # For initial frame, normalize by search_radius
                end
                if norm_dist_factor < 1f-5 norm_dist_factor = 1f0 end # Avoid division by zero

                for (k, cand_feat) in enumerate(candidate_features_this_frame)
                    dist_to_search_center = distance_2d(cand_feat.centroid_px, search_center)
                    
                    # Normalized distance score (0 to ~1, can be >1 if outside WORM_DISTANCE_LINK but still a candidate by other means)
                    dist_score = dist_to_search_center / norm_dist_factor 
                    
                    # Circularity score (0 to 1, where 0 is perfect circularity of 1.0)
                    circ_score = 1.0f0 - cand_feat.circularity
                    
                    # Combined score (example weights: 60% distance, 40% circularity)
                    # Adjust weights as needed.
                    current_score = 0.6f0 * dist_score + 0.4f0 * circ_score
                    
                    if current_score < best_score
                        best_score = current_score
                        best_candidate_idx = k
                    end
                end
                if best_candidate_idx != -1
                    found_feature_this_frame = candidate_features_this_frame[best_candidate_idx]
                end
            end
            
            # Final validity checks (arena boundary and absolute circularity)
            if !isnothing(found_feature_this_frame) &&
               distance_2d(found_feature_this_frame.centroid_px, img_center_px) > arena_radius_px + arena_boundary_tolerance_px
                found_feature_this_frame = nothing 
            end
            if !isnothing(found_feature_this_frame) &&
               found_feature_this_frame.circularity < MIN_ACCEPTABLE_CIRCULARITY
                found_feature_this_frame = nothing 
            end
        end

        pass_features_output[current_frame_idx] = found_feature_this_frame
        if !isnothing(found_feature_this_frame)
            _last_identified_feature_this_pass = found_feature_this_frame
            _last_pos_for_search_expansion_this_pass = found_feature_this_frame.centroid_px
            _consecutive_misses_this_pass = 0
        else
            _consecutive_misses_this_pass += 1
            _last_identified_feature_this_pass = nothing
        end
    end
    return pass_features_output
end

# ... (_core_worm_tracking_bidirectional_final_interpolate remains the same as it calls the modified single_pass)
function _core_worm_tracking_bidirectional_final_interpolate(
    binary_stack_global::BitArray{3},
    num_frames::Int,
    img_rows::Int, img_cols::Int
)
    img_center_px = Point2f(img_cols / 2.0f0, img_rows / 2.0f0)
    arena_radius_px = min(img_rows, img_cols) / 2.1f0
    on_circle_tolerance_px = 15.0f0; arena_boundary_tolerance_px = 10.0f0 # Arena boundary for final check
    
    # Call the (now modified) single pass function
    forward_pass_features = _core_worm_tracking_single_pass(binary_stack_global, num_frames, img_rows, img_cols, true, img_center_px, arena_radius_px, on_circle_tolerance_px, arena_boundary_tolerance_px)
    backward_pass_features = _core_worm_tracking_single_pass(binary_stack_global, num_frames, img_rows, img_cols, false, img_center_px, arena_radius_px, on_circle_tolerance_px, arena_boundary_tolerance_px)
    
    merged_features = Vector{Union{Nothing, WormFeatures}}(nothing, num_frames)
    for i in 1:num_frames
        f_feat = forward_pass_features[i]; b_feat = backward_pass_features[i]
        if !isnothing(f_feat) merged_features[i] = f_feat
        elseif !isnothing(b_feat) merged_features[i] = b_feat
        else merged_features[i] = nothing end
    end
    positions_px_out = Point3f[]; areas_px2_out = Int[]; perimeters_px_out = Float32[]; circularities_out = Float32[]
    frame_idx = 1
    while frame_idx <= num_frames
        current_feature = merged_features[frame_idx]
        if !isnothing(current_feature)
            push!(positions_px_out, Point3f(current_feature.centroid_px[1], current_feature.centroid_px[2], Float32(frame_idx)))
            push!(areas_px2_out, current_feature.area_px2); push!(perimeters_px_out, current_feature.perimeter_px); push!(circularities_out, current_feature.circularity)
            frame_idx += 1
        else
            gap_start_frame = frame_idx
            while frame_idx <= num_frames && isnothing(merged_features[frame_idx]) frame_idx += 1 end
            feature_before_gap = gap_start_frame > 1 ? merged_features[gap_start_frame - 1] : nothing
            feature_after_gap = frame_idx <= num_frames ? merged_features[frame_idx] : nothing
            if !isnothing(feature_before_gap) && !isnothing(feature_after_gap)
                num_intermediate_points = (frame_idx - 1) - gap_start_frame + 1
                if num_intermediate_points > 0
                    interpolated_points_2d = interpolate_circular_path(feature_before_gap.centroid_px, feature_after_gap.centroid_px, img_center_px, arena_radius_px, num_intermediate_points)
                    for k in 1:num_intermediate_points
                        interp_frame_actual_idx = gap_start_frame + k - 1
                        if k <= length(interpolated_points_2d)
                            interp_pos_2d = interpolated_points_2d[k]
                            push!(positions_px_out, Point3f(interp_pos_2d[1], interp_pos_2d[2], Float32(interp_frame_actual_idx)))
                            push!(areas_px2_out, feature_before_gap.area_px2); push!(perimeters_px_out, feature_before_gap.perimeter_px); push!(circularities_out, feature_before_gap.circularity)
                        end
                    end
                end
            end
        end
    end
    return positions_px_out, areas_px2_out, perimeters_px_out, circularities_out, img_rows, img_cols
end

# --- Plotting Functions (Unchanged) ---
# ... (plot_all_distances_from_center, calculate_msd, plot_all_msd, plot_sum_pairwise_differences_over_time)
# ... (view_binary_stack_with_tracks_glmakie)
function plot_all_distances_from_center(all_worm_data::Vector{WormProcessedData})
    fig = Figure(size=(1000, 750))
    ax = Axis(fig[1,1], xlabel="Time (s)", ylabel="Distance from Center (cm)", title="Distance from Image Center Over Time")
    colors = distinguishable_colors(max(1,length(all_worm_data)), [RGB(1,1,1), RGB(0,0,0)], lchoices=range(20, stop=70, length=15))
    for (i, data) in enumerate(all_worm_data)
        center_x_px = data.img_cols / 2.0f0; center_y_px = data.img_rows / 2.0f0
        center_point_px = Point2f(center_x_px, center_y_px)
        distances_cm = Float32[]
        valid_times_s = data.times_s
        for p_idx in eachindex(data.positions_px)
            p_px = data.positions_px[p_idx]
            current_pos_2d_px = Point2f(p_px[1], p_px[2])
            dist_px = distance_2d(current_pos_2d_px, center_point_px)
            push!(distances_cm, dist_px * CM_PER_PIXEL)
        end
        lines!(ax, valid_times_s, distances_cm, color=colors[i], linewidth=2, label=data.name)
    end
    if length(all_worm_data) > 0 axislegend(ax) end
    display(GLMakie.Screen(), fig); return fig
end

function calculate_msd(positions_2d_cm::Vector{Point2f}, times_s::Vector{Float32})
    n_points = length(positions_2d_cm)
    if n_points < 2 || length(times_s) != n_points return Dict{Float32, Float32}() end
    msd_values = Dict{Float32, Vector{Float32}}()
    for lag_frames in 1:n_points-1
        for i in 1:(n_points - lag_frames)
            p1 = positions_2d_cm[i]; p2 = positions_2d_cm[i + lag_frames]
            diff_vec = p2 - p1; dist_sq = sum(diff_vec.^2)
            time_lag_s = abs(times_s[i + lag_frames] - times_s[i])
            time_lag_key = round(time_lag_s, digits=3)
            if !haskey(msd_values, time_lag_key) msd_values[time_lag_key] = Float32[] end
            push!(msd_values[time_lag_key], dist_sq)
        end
    end
    averaged_msd = Dict{Float32, Float32}()
    for (lag, sds) in msd_values
        if !isempty(sds) averaged_msd[lag] = mean(sds) end
    end
    return averaged_msd
end

function plot_all_msd(all_worm_data::Vector{WormProcessedData})
    fig = Figure(size=(1000, 750))
    ax = Axis(fig[1,1], xlabel="Time Lag τ (s)", ylabel="MSD (cm²)", xscale=log10, yscale=log10)
    colors = distinguishable_colors(max(1,length(all_worm_data)), [RGB(1,1,1), RGB(0,0,0)], lchoices=range(20, stop=70, length=15))
    min_msd_val_overall = Inf32; max_msd_val_overall = -Inf32; min_lag_overall = Inf32; max_lag_overall = -Inf32
    plotted_anything = false
    for (i, data) in enumerate(all_worm_data)
        if isempty(data.positions_cm) || isempty(data.times_s) || length(data.positions_cm) != length(data.times_s) continue end
        positions_2d_cm = [Point2f(p[1], p[2]) for p in data.positions_cm]
        msd_dict = calculate_msd(positions_2d_cm, data.times_s)
        if isempty(msd_dict) continue end
        sorted_lags_s = sort(collect(keys(msd_dict))); sorted_msd_values = [msd_dict[lag] for lag in sorted_lags_s]
        valid_indices = (sorted_lags_s .> 1e-9) .& (sorted_msd_values .> 1e-9)
        plot_lags = sorted_lags_s[valid_indices]; plot_msd = sorted_msd_values[valid_indices]
        if !isempty(plot_lags) && !isempty(plot_msd)
            lines!(ax, plot_lags, plot_msd, color=colors[i], linewidth=2, label=data.name)
            min_msd_val_overall = min(min_msd_val_overall, minimum(plot_msd)); max_msd_val_overall = max(max_msd_val_overall, maximum(plot_msd))
            min_lag_overall = min(min_lag_overall, minimum(plot_lags)); max_lag_overall = max(max_lag_overall, maximum(plot_lags))
            plotted_anything = true
        end
    end
    if plotted_anything && isfinite(min_lag_overall) && isfinite(max_lag_overall) && min_lag_overall < max_lag_overall
        actual_min_lag = max(min_lag_overall, 1e-3f0); actual_max_lag = max_lag_overall
        if actual_min_lag < actual_max_lag
            ref_lags = 10 .^ LinRange(log10(actual_min_lag), log10(actual_max_lag), 50)
            anchor_y_at_min_lag = isfinite(min_msd_val_overall) ? max(min_msd_val_overall, 1e-9f0) : 1e-4f0
            y1_vals = anchor_y_at_min_lag .* (ref_lags ./ ref_lags[1]).^1.0; lines!(ax, ref_lags, y1_vals, color=RGB(0.5,0.5,0.5), linestyle=:dash, label="α=1")
            y2_vals = anchor_y_at_min_lag .* (ref_lags ./ ref_lags[1]).^2.0; lines!(ax, ref_lags, y2_vals, color=RGB(0.7,0.7,0.7), linestyle=:dot, label="α=2")
            y05_vals = anchor_y_at_min_lag .* (ref_lags ./ ref_lags[1]).^0.5; lines!(ax, ref_lags, y05_vals, color=RGB(0.3,0.3,0.3), linestyle=:dashdot, label="α=0.5")
        end
    end
    if plotted_anything axislegend(ax, position=:rb) end
    display(GLMakie.Screen(), fig); return fig
end

function plot_sum_pairwise_differences_over_time(all_worm_data::Vector{WormProcessedData})
    fig = Figure(size=(1000, 750))
    ax = Axis(fig[1,1], xlabel="Time (s)", ylabel="Sum of Pairwise Distances (cm)", title="Sum of Pairwise Worm Distances Over Time")
    if length(all_worm_data) < 2 display(GLMakie.Screen(), fig); return fig end
    all_frame_indices_with_data = Set{Int}(); frame_to_time_map = Dict{Int, Float32}()
    worm_pos_at_frame_idx = [Dict{Int, Point2f}() for _ in 1:length(all_worm_data)]
    for (data_idx, data) in enumerate(all_worm_data)
        for p3d_px in data.positions_px
            frame_idx = round(Int, p3d_px[3]); time_s = frame_idx * MSEC_PER_FRAME / 1000.0f0
            push!(all_frame_indices_with_data, frame_idx)
            if !haskey(frame_to_time_map, frame_idx) frame_to_time_map[frame_idx] = time_s end
            pos_x_cm = p3d_px[1] * CM_PER_PIXEL; pos_y_cm = (data.img_rows - p3d_px[2]) * CM_PER_PIXEL
            worm_pos_at_frame_idx[data_idx][frame_idx] = Point2f(pos_x_cm, pos_y_cm)
        end
    end
    if isempty(all_frame_indices_with_data) display(GLMakie.Screen(), fig); return fig end
    sorted_common_frame_indices = sort(collect(all_frame_indices_with_data))
    sum_diff_values = Float32[]; valid_times_for_sum_diff = Float32[]
    for current_frame_idx in sorted_common_frame_indices
        active_worms_positions_at_frame = Point2f[]
        for data_idx in 1:length(all_worm_data)
            if haskey(worm_pos_at_frame_idx[data_idx], current_frame_idx)
                push!(active_worms_positions_at_frame, worm_pos_at_frame_idx[data_idx][current_frame_idx])
            end
        end
        if length(active_worms_positions_at_frame) >= 2
            current_sum_pairwise_dist = 0.0f0
            for i in 1:length(active_worms_positions_at_frame), j in (i+1):length(active_worms_positions_at_frame)
                current_sum_pairwise_dist += distance_2d(active_worms_positions_at_frame[i], active_worms_positions_at_frame[j])
            end
            push!(sum_diff_values, current_sum_pairwise_dist)
            push!(valid_times_for_sum_diff, frame_to_time_map[current_frame_idx])
        end
    end
    lines!(ax, valid_times_for_sum_diff, sum_diff_values, color=:green, linewidth=2)
    display(GLMakie.Screen(), fig); return fig
end

function view_binary_stack_with_tracks_glmakie(
    binary_stack_to_show::BitArray{3},
    tracked_worm_data::WormProcessedData;
    marker_color=RGBAf(1.0, 0.0, 0.0, 0.7) 
)
    rows, cols, num_total_frames = size(binary_stack_to_show)
    worm_info_for_frame = Dict{Int, Tuple{Point2f, Int}}() 
    for i in 1:length(tracked_worm_data.positions_px)
        pt3 = tracked_worm_data.positions_px[i]
        frame_idx = Int(round(pt3[3]))
        centroid = Point2f(pt3[2], pt3[1]) # (col_px, row_px) for GLMakie point
        area = tracked_worm_data.areas_px2[i]
        worm_info_for_frame[frame_idx] = (centroid, area)
    end
    fig = Figure(size = (cols > rows ? (800, 800 * rows/cols + 100) : (800 * cols/rows + 100, 800)))
    ax_img = Axis(fig[1, 1])
    slider = Slider(fig[2, 1], range=1:num_total_frames, startvalue=1)
    frame_idx_obs = slider.value
    img_slice_obs = lift(frame_idx_obs) do f_idx
        return Gray.(view(binary_stack_to_show, :, :, f_idx)) 
    end
    image!(ax_img, img_slice_obs, interpolate=false, colormap=:grays, colorrange=(0,1))
    worm_marker_obs = Observable([Circle(Point2f(NaN, NaN), 0f0)])
    on(frame_idx_obs) do f_idx
        if haskey(worm_info_for_frame, f_idx)
            centroid_px, area_px2 = worm_info_for_frame[f_idx]
            radius = sqrt(max(0.0, Float64(area_px2)) / π)
            display_radius = clamp(Float32(radius), 2f0, 25f0)
            worm_marker_obs[] = [Circle(centroid_px, display_radius)]
        else
            worm_marker_obs[] = [Circle(Point2f(NaN, NaN), 0f0)]
        end
    end
    poly!(ax_img, worm_marker_obs, color=marker_color, strokecolor=:transparent)
    display(GLMakie.Screen(), fig)
    return fig
end


filepaths = [
     "/Users/dl0346/Documents/PlanarianVideos/May28/TrailFollow/ATM1_1.tif",
          "/Users/dl0346/Documents/PlanarianVideos/May28/TrailFollow/ATM1_2.tif",
               "/Users/dl0346/Documents/PlanarianVideos/May28/TrailFollow/ATM1_3.tif",
]
global binary_stack_for_debug = BitArray{3}(undef, (0,0,0))
all_processed_data = WormProcessedData[]

for (idx, filepath) in enumerate(filepaths) 
    println("Processing: $filepath")
    raw_stack_data = TiffImages.load(filepath)
    img_rows, img_cols = size(raw_stack_data, 1), size(raw_stack_data, 2)
    num_raw_frames = size(raw_stack_data, 3)
    float_stack_inverted = 1.0f0 .- Float32.(raw_stack_data)
    stabilized_stack_for_processing = float_stack_inverted
    stabilize_flag = false
    if stabilize_flag
        reference_frame_stabilization = median(float_stack_inverted, dims=3)[:,:,1]
        stabilized_stack_temp = similar(float_stack_inverted)
        frames_to_process = 1:num_raw_frames
        @showprogress Threads.@threads for i in frames_to_process
            current_frame = float_stack_inverted[:, :, i]
            dy, dx = phase_correlate(reference_frame_stabilization, current_frame)
            stabilized_stack_temp[:, :, i] = translate_image(current_frame, -dy, -dx, 0.0f0)
        end
        stabilized_stack_for_processing = stabilized_stack_temp
    end
    new_ref_division = median(stabilized_stack_for_processing, dims=3)
    final_processed_stack = stabilized_stack_for_processing ./ (new_ref_division .+ eps(Float32(1.0)))
    global binary_stack_for_debug = generate_full_binary_stack(final_processed_stack, THRESHOLD_VALUE)
    final_num_frames_to_track = size(binary_stack_for_debug, 3)
    positions_px, areas_px2, perimeters_px, circularities_data, track_img_rows, track_img_cols =
        _core_worm_tracking_bidirectional_final_interpolate(binary_stack_for_debug, final_num_frames_to_track, img_rows, img_cols)
    video_short_name = basename(filepath)
    times_s_vec = [p[3] * MSEC_PER_FRAME / 1000.0f0 for p in positions_px]
    positions_cm_s_vec = [Point3f(p[1]*CM_PER_PIXEL, (track_img_rows - p[2])*CM_PER_PIXEL, times_s_vec[i]) for (i,p) in enumerate(positions_px)]
    data = WormProcessedData(filepath, video_short_name, positions_px, areas_px2, perimeters_px, circularities_data, img_rows, img_cols, times_s_vec, positions_cm_s_vec)
    push!(all_processed_data, data)
end

## Original GLMakie 3D Plotting Section (Trajectory Plot)
all_worm_data_for_glmakie = all_processed_data
fig_3d_traj = Figure(size=(1200, 900)) 
max_time_s_overall = 0.0f0; img_width_cm_overall = 7.5f0; img_height_cm_overall = 7.5f0
valid_times_all_data = [data.times_s for data in all_worm_data_for_glmakie if !isempty(data.times_s)]

max_times_per_dataset = [maximum(ts) for ts in valid_times_all_data if !isempty(ts)]
first_data_for_dims = all_worm_data_for_glmakie[1]
img_width_cm_overall = first_data_for_dims.img_cols * CM_PER_PIXEL
img_height_cm_overall = first_data_for_dims.img_rows * CM_PER_PIXEL

ax_3d = Axis3(fig_3d_traj[1,1], xlabel="X (cm)", ylabel="Y (cm)", zlabel="Time (s)", 
limits = ( (0, img_width_cm_overall),(0, img_height_cm_overall),(0, maximum(max_times_per_dataset))  ))

colors_traj = distinguishable_colors(max(1,length(all_worm_data_for_glmakie)), [RGB(0,0,1), RGB(0,0,0)], lchoices=range(20, stop=70, length=15))
for (i, data) in enumerate(all_worm_data_for_glmakie)
    lines!(ax_3d, data.positions_cm, color=colors_traj[i], linewidth=4, label=data.name)
end
display(GLMakie.Screen(), fig_3d_traj) 

#data_for_visualization = all_processed_data[end] 
#interactive_fig = view_binary_stack_with_tracks_glmakie(binary_stack_for_debug, data_for_visualization)

