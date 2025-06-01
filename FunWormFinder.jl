module FUNC_WormFinder

using LinearAlgebra # For norm (distance calculation)
using Statistics    # For mean
using GeometryBasics # For Polygon, Rect, Point2f, Point3f, Circle\
include("FUNC_ObjectFind.jl")
using .FUNC_ObjectFind

export find_objects_features, create_track, worm_track_pass, WormData

# Constants for worm tracking
const MIN_AREA = 15
const WORM_DISTANCE_LINK = 40.0f0
const WORM_DISTANCE_SEARCHADD = 0.1f0
const MIN_ACCEPTABLE_CIRCULARITY = 0.10f0
const MAX_AREA = 500

# Helper distance function
distance_2d(p1::Point2f, p2::Point2f) = norm(p1 - p2)

# Features for each detected worm
struct WormFeatures
    centroid::Point2f
    wormarea::Int
    major_axis::Float32
    minor_axis::Float32
    pixels_abs::Vector{Tuple{Int,Int}} # (row, col)
    boundary_pixels_abs::Vector{Point2f}
end

# Data container for full worm tracking results
struct WormData
    filepath::String
    positions_::Vector{Point3f}
    areas_2::Vector{Int}
    major_axes::Vector{Float32}
    minor_axes::Vector{Float32}
    img_rows::Int
    img_cols::Int
    times_s::Vector{Float32}
    positions_cm::Vector{Point3f}
end

# Single pass worm tracking
function worm_track_pass(binary, is_forward_pass, img_center, arena_radius, on_circle_tolerance)::Vector{Union{Nothing, WormFeatures}}
    img_rows, img_cols, num_frames = size(binary, 1), size(binary, 2), size(binary, 3)
    worm_track = Vector{Union{Nothing, WormFeatures}}(nothing, num_frames)
    last_feature::Union{Nothing, WormFeatures} = nothing
    last_position::Union{Nothing, Point2f} = nothing
    num_misses = 0

    frame_iterator = is_forward_pass ? (1:num_frames) : (num_frames:-1:1)
    for current_frame in frame_iterator
        chosen_worm::Union{Nothing, WormFeatures} = nothing
        possible_worm = WormFeatures[]
        search_center::Point2f = Point2f(0,0)

        binary_view = view(binary, :, :, current_frame)

        
        all_objects_in_frame = [WormFeatures(f.centroid, f.wormarea, f.major_axis, f.minor_axis, f.pixels_abs, f.boundary_pixels_abs) for f in find_objects_features(binary_view, MIN_AREA, MAX_AREA)]

        if !isnothing(last_feature)
            search_center = last_feature.centroid
            search_radius = WORM_DISTANCE_LINK
            possible_worm = filter(f -> distance_2d(f.centroid, search_center) <= search_radius, all_objects_in_frame)
        elseif !isnothing(last_position)
            search_center = last_position
            search_radius = WORM_DISTANCE_LINK + num_misses * WORM_DISTANCE_SEARCHADD
            potential_candidates = filter(f -> distance_2d(f.centroid, search_center) <= search_radius, all_objects_in_frame)
            possible_worm = filter(f -> abs(distance_2d(f.centroid, img_center) - arena_radius) <= on_circle_tolerance, potential_candidates)
        else
            search_center = img_center
            search_radius = arena_radius
            possible_worm = filter(f -> distance_2d(f.centroid, img_center) <= search_radius, all_objects_in_frame)
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
            if !isnothing(chosen_worm) && chosen_worm.circularity < MIN_ACCEPTABLE_CIRCULARITY
                chosen_worm = nothing
            end
        end

        worm_track[current_frame] = chosen_worm
        if !isnothing(chosen_worm)
            last_feature = chosen_worm
            last_position = chosen_worm.centroid
            num_misses = 0
        else
            num_misses += 1
            last_feature = nothing
        end
    end
    return worm_track
end

# Create full track with interpolation for gaps
function create_track(binary)
    img_rows, img_cols, num_frames = size(binary, 1), size(binary, 2), size(binary, 3)
    img_center = Point2f(img_cols / 2.0f0, img_rows / 2.0f0)
    arena_radius = min(img_rows, img_cols) / 2.1f0
    on_circle_tolerance = 15.0f0

    forward_pass_features = worm_track_pass(binary, true, img_center, arena_radius, on_circle_tolerance)
    backward_pass_features = worm_track_pass(binary, false, img_center, arena_radius, on_circle_tolerance)

    merged_features = Vector{Union{Nothing, WormFeatures}}(nothing, num_frames)
    for i in 1:num_frames
        f_feat = forward_pass_features[i]
        b_feat = backward_pass_features[i]
        if !isnothing(f_feat)
            merged_features[i] = f_feat
        elseif !isnothing(b_feat)
            merged_features[i] = b_feat
        else
            merged_features[i] = nothing
        end
    end

    positions = Point3f[]
    areas = Int[]
    major_axes = Float32[]
    minor_axes = Float32[]
    frame_idx = 1

    while frame_idx <= num_frames
        worm = merged_features[frame_idx]
        if !isnothing(worm)
            push!(positions, Point3f(worm.centroid[1], worm.centroid[2], Float32(frame_idx)))
            push!(areas, worm.wormarea)
            push!(major_axes, worm.major_axis)
            push!(minor_axes, worm.minor_axis)
            frame_idx += 1
        else
            gap_start_frame = frame_idx
            while frame_idx <= num_frames && isnothing(merged_features[frame_idx])
                frame_idx += 1
            end
            feature_before_gap = gap_start_frame > 1 ? merged_features[gap_start_frame - 1] : nothing
            feature_after_gap = frame_idx <= num_frames ? merged_features[frame_idx] : nothing
            if !isnothing(feature_before_gap) && !isnothing(feature_after_gap)
                num_intermediate_points = frame_idx - gap_start_frame
                interpolated_points_2d = interpolate_circular_path(feature_before_gap.centroid, feature_after_gap.centroid, img_center, arena_radius, num_intermediate_points)
                for k in 1:num_intermediate_points
                    interp_frame_actual_idx = gap_start_frame + k - 1
                    interp_pos_2d = interpolated_points_2d[k]
                    push!(positions, Point3f(interp_pos_2d[1], interp_pos_2d[2], Float32(interp_frame_actual_idx)))
                    push!(areas, feature_before_gap.wormarea)
                    push!(major_axes, feature_before_gap.major_axis)
                    push!(minor_axes, feature_before_gap.minor_axis)
                end
            end
        end
    end
    return positions, areas, major_axes, minor_axes, img_rows, img_cols
end

end # module FUNC_WormFinderCIRCLE