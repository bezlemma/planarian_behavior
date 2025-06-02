#When doing a circle assay, we have special information that could help us find the worm at the edges of the circle
#TODO: Make this pass variables rather than using constants
#TODO: Make some logic about primary axis being most likely axis that worm moves in

module FUNC_WormFinder

using LinearAlgebra, Statistics, GeometryBasics 

include("FUNC_ObjectFind.jl"); using .FUNC_ObjectFind

export find_objects_features, create_track, worm_track_pass, interpolate_circular_path, WormData

# Constants for worm tracking
const WORM_DISTANCE_LINK = 40.0

# Helper distance function
distance_2d(p1::Point2f, p2::Point2f) = norm(p1 - p2)

# Features for each detected worm
struct WormFeatures
    centroid::Point2f
    wormarea::Int
    major_axis::Float32; minor_axis::Float32
    pixels_abs::Vector{Tuple{Int,Int}} # (row, col)
    boundary_pixels_abs::Vector{Point2f}
end

# Data container for full worm tracking results
struct WormData
    filepath::String
    positions_::Vector{Point3f}
    areas_2::Vector{Int}
    major_axes::Vector{Float32}; minor_axes::Vector{Float32}
    img_rows::Int; img_cols::Int
    times_s::Vector{Float32}; positions_cm::Vector{Point3f}
    behaviors::Vector{Symbol}  # new field for behaviors
end

# Single pass worm tracking
function worm_track_pass(binary, is_forward_pass, img_center, arena_radius, on_circle_tol, WORM_DISTANCE_SEARCHADD)
    
    num_frames = size(binary, 3)
    worm_track = Vector{Union{Nothing, WormFeatures}}(nothing, num_frames)
    last_feature::Union{Nothing, WormFeatures} = nothing
    last_position::Union{Nothing, Point2f} = nothing
    num_misses = 0

    frame_iterator = is_forward_pass ? (1:num_frames) : (num_frames:-1:1)
    for current_frame in frame_iterator
        chosen_worm::Union{Nothing, WormFeatures} = nothing
        possible_worm = WormFeatures[]
        search_center = Point2f(0,0)

        binary_view = view(binary, :, :, current_frame)
        all_objects = [WormFeatures(f.centroid, f.wormarea, f.major_axis, f.minor_axis, f.pixels_abs, f.boundary_pixels_abs) for f in find_objects_features(binary_view)]

        if last_feature != nothing #Do we have a feature from the last frame
 
            possible_worm = filter(f -> 
            distance_2d(f.centroid, last_feature.centroid) <= WORM_DISTANCE_LINK, all_objects)

        elseif last_position != nothing #If not, work off the last known position
            search_radius = WORM_DISTANCE_LINK + num_misses * WORM_DISTANCE_SEARCHADD

            potential_candidates = filter(f -> 
            distance_2d(f.centroid, last_position) <= search_radius, all_objects)

            possible_worm = filter(f -> 
            abs(distance_2d(f.centroid, img_center) - arena_radius) <= on_circle_tol, potential_candidates)
        
        else #if we have nothing, pick the worm candidate closest to the center.
            search_center = img_center
            search_radius = arena_radius
            possible_worm = filter(f -> 
            distance_2d(f.centroid, img_center) <= search_radius, all_objects)
        end

        if !isempty(possible_worm)
            if length(possible_worm) == 1
                chosen_worm = possible_worm[1]
            else
                best_score = Inf32
                best_candidate_idx = -1
                norm_dist_factor = WORM_DISTANCE_LINK
                for (k, cand_feat) in enumerate(possible_worm)
                    dist_to_search_center = distance_2d(cand_feat.centroid, search_center)
                    dist_score = dist_to_search_center / norm_dist_factor
                    # score based solely on distance
                    current_score = dist_to_search_center / norm_dist_factor
                    if current_score < best_score
                        best_score = current_score
                        best_candidate_idx = k
                    end
                end
                if best_candidate_idx != -1
                    chosen_worm = possible_worm[best_candidate_idx]
                end
            end
            worm_track[current_frame] = chosen_worm

            #Save our worm or lack of worm for next steps search logic
            if chosen_worm != nothing
                last_feature = chosen_worm
                last_position = chosen_worm.centroid
                num_misses = 0
            else
                num_misses += 1
                last_feature = nothing
            end

        end
    end
    return worm_track
end

# Create full track with interpolation for gaps
function create_track(binary, circle, WORM_DISTANCE_SEARCHADD)

    #Set up image info
    img_rows, img_cols, num_frames = size(binary, 1), size(binary, 2), size(binary, 3)
    img_center = Point2f(img_cols / 2.0f0, img_rows / 2.0f0)
    arena_radius = min(img_rows, img_cols) / 2
    on_circle_tol = 15

    #Look for the worm both forwards and backwards
    forward_pass_features = worm_track_pass(binary, true, img_center, arena_radius, on_circle_tol,WORM_DISTANCE_SEARCHADD)
    backward_pass_features = worm_track_pass(binary, false, img_center, arena_radius, on_circle_tol,WORM_DISTANCE_SEARCHADD)

    merged_features = Vector{Union{Nothing, WormFeatures}}(nothing, num_frames)
    positions = Point3f[]
    areas = Int[]
    major_axes = Float32[]; minor_axes = Float32[];
    frame_idx = 1

    #Use the forward, but if there is no forward use the backwards
    for i in 1:num_frames
        f_feat = forward_pass_features[i]
        b_feat = backward_pass_features[i]
        if f_feat != nothing
            merged_features[i] = f_feat
        elseif b_feat != nothing
            merged_features[i] = b_feat
        else
            merged_features[i] = nothing
        end
    end

    #
    
    while frame_idx <= num_frames
        worm = merged_features[frame_idx]
        if worm != nothing #Ok, we have a worm, just add it.
            push!(positions, Point3f(worm.centroid[1], worm.centroid[2], Float32(frame_idx)))
            push!(areas, worm.wormarea)
            push!(major_axes, worm.major_axis); push!(minor_axes, worm.minor_axis)
            frame_idx += 1
        elseif circle && worm == nothing #If we know the worm exists in a circle, we can guess where it is on the circle     
            gap_start_frame = frame_idx
            
            #find the bounds of the gap
            while (frame_idx <= num_frames) && (merged_features[frame_idx] == nothing)
                frame_idx += 1
            end

            #Grab the worm on each side of the gap
            worm_before_gap = gap_start_frame > 1 ? merged_features[gap_start_frame - 1] : nothing
            worm_after_gap = frame_idx <= num_frames ? merged_features[frame_idx] : nothing

            #It's possible that there is no worm at the end of the gap, but if there is, we can interpolate
            if worm_before_gap != nothing && worm_after_gap != nothing

                num_intermediate_points = frame_idx - gap_start_frame
                
                #Use the mean radius of the worm position from the center of the image before and after the gap
                dist_before_gap = distance_2d(worm_before_gap.centroid, img_center)
                dist_after_gap = distance_2d(worm_after_gap.centroid, img_center)
                radius = (dist_before_gap + dist_after_gap) / 2

                interpolated_points_2d = interpolate_circular_path(
                    worm_before_gap.centroid,
                    worm_after_gap.centroid,
                    img_center,
                    radius,
                    num_intermediate_points,
                    WORM_DISTANCE_SEARCHADD)
                # only proceed if interpolation returned the full set of points
                if length(interpolated_points_2d) == num_intermediate_points
                    for k in 1:num_intermediate_points
                        interp_frame_actual_idx = gap_start_frame + k - 1
                        interp_pos_2d = interpolated_points_2d[k]
                        push!(positions, Point3f(interp_pos_2d[1], interp_pos_2d[2], Float32(interp_frame_actual_idx)))
                        push!(areas, worm_before_gap.wormarea)
                        push!(major_axes, worm_before_gap.major_axis)
                        push!(minor_axes, worm_before_gap.minor_axis)
                    end
                end

            end

        end
    end

    #  post‐filter long jumps ---
    begin
        filtered_positions = Point3f[]
        filtered_areas     = Int[]
        filtered_major_axes = Float32[]
        filtered_minor_axes = Float32[]
        if !isempty(positions)
            # always keep first
            push!(filtered_positions, positions[1]); push!(filtered_areas, areas[1])
            push!(filtered_major_axes, major_axes[1]); push!(filtered_minor_axes, minor_axes[1])
            for i in 2:length(positions)
                prev = filtered_positions[end]
                curr = positions[i]
                if distance_2d(Point2f(prev[1], prev[2]), Point2f(curr[1], curr[2])) <= WORM_DISTANCE_SEARCHADD
                    push!(filtered_positions, curr)
                    push!(filtered_areas, areas[i])
                    push!(filtered_major_axes, major_axes[i])
                    push!(filtered_minor_axes, minor_axes[i])
                end
            end
            positions    = filtered_positions
            areas         = filtered_areas
            major_axes   = filtered_major_axes
            minor_axes   = filtered_minor_axes
        end
    end

    return positions, areas, major_axes, minor_axes
end

# Interpolate points along circular arc between two points on circle
function interpolate_circular_path(p_start_, p_end_, circle_center_, circle_radius_, num_intermediate_points, WORM_DISTANCE_SEARCHADD)::Vector{Point2f}
    v_start = p_start_ - circle_center_; v_end = p_end_ - circle_center_
    angle_start_rad = atan(v_start[2], v_start[1]); angle_end_rad = atan(v_end[2], v_end[1])
    delta_angle_rad = angle_end_rad - angle_start_rad
    if delta_angle_rad > π delta_angle_rad -= 2π elseif delta_angle_rad < -π delta_angle_rad += 2π end
    interpolated_points = Vector{Point2f}(undef, num_intermediate_points)
    total_steps_for_interpolation = num_intermediate_points + 1
    for i in 1:num_intermediate_points
        fraction = Float32(i) / Float32(total_steps_for_interpolation)
        current_angle_rad = angle_start_rad + fraction * delta_angle_rad
        pt_x = circle_center_[1] + circle_radius_ * cos(current_angle_rad)
        pt_y = circle_center_[2] + circle_radius_ * sin(current_angle_rad)
        interpolated_points[i] = Point2f(pt_x, pt_y)
    end

    # reject the entire path if any per‐step jump exceeds allowed threshold
    prev_pt = p_start_
    for pt in interpolated_points
        if distance_2d(prev_pt, pt) > WORM_DISTANCE_SEARCHADD
            return Point2f[]
        end
        prev_pt = pt
    end
    if distance_2d(prev_pt, p_end_) > WORM_DISTANCE_SEARCHADD
        return Point2f[]
    end

    return interpolated_points
end

end # module FUNC_WormFinder