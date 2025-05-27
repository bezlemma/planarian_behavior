using TiffImages
using GLMakie
using LinearAlgebra # For norm (distance calculation)
using Statistics    # For mean
using Colors        # For distinguishable_colors, RGB, Gray
using GeometryBasics # For Polygon

# --- Configuration ---
const TIFF_FILE_PATH = "/Users/dl0346/Documents/PlanarianVideos/ColloidTestInverse.tif"
const THRESHOLD_VALUE = 0.5      # Grayscale threshold (0.0 to 1.0)
const MIN_WORM_AREA = 10         # Minimum pixel area to be considered a worm
const MAX_LINKING_DISTANCE = 100.0 # Max pixels a worm can move between frames to be linked

# --- Physical Unit Conversion ---
const PIXELS_PER_CM = 58.6 # pixels/cm
const MSEC_PER_FRAME = 50.0 * 2.0 # msec per frame

# --- Helper Functions ---
distance_2d(p1, p2) = norm(p1 - p2)

struct WormFeatures
    centroid::Point2f
    area::Int
    perimeter::Float32
    circularity::Float32
    pixels::Vector{Tuple{Int,Int}} # Absolute pixel coordinates (r, c)
    boundary_pixels::Vector{Point2f} # For plotting outline (x, y)
end

"""
Finds connected components (objects) in a binary mask using BFS.
Calculates centroid, area, perimeter, and circularity for each object.
"""
function find_objects_features(binary_mask::BitMatrix, min_area::Int)
    rows, cols = size(binary_mask)
    visited = falses(rows, cols)
    object_features_list = WormFeatures[]

    NEIGHBORS_8 = ((-1, -1), (-1, 0), (-1, 1), (0, -1), (0, 1),(1, -1), (1, 0), (1, 1))
    NEIGHBORS_4_FOR_PERIMETER = ((-1,0), (1,0), (0,-1), (0,1))

    for r_start in 1:rows
        for c_start in 1:cols
            if binary_mask[r_start, c_start] && !visited[r_start, c_start]
                current_object_pixels_tuples = Tuple{Int,Int}[]
                q = Tuple{Int,Int}[(r_start, c_start)]
                visited[r_start, c_start] = true
                head = 1

                while head <= length(q)
                    curr_r, curr_c = q[head]
                    head += 1
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

                obj_area = length(current_object_pixels_tuples)
                if obj_area >= min_area
                    sum_x, sum_y = 0.0f0, 0.0f0
                    for (pix_r, pix_c) in current_object_pixels_tuples
                        sum_x += pix_c
                        sum_y += pix_r
                    end
                    centroid_x = sum_x / obj_area
                    centroid_y = sum_y / obj_area

                    # Calculate perimeter (count of exposed edges) and boundary pixels
                    obj_perimeter_edges = 0
                    boundary_pixel_coords = Point2f[]
                    pixel_set = Set(current_object_pixels_tuples)

                    for (pix_r, pix_c) in current_object_pixels_tuples
                        is_on_boundary_for_plot = false
                        for (dr, dc) in NEIGHBORS_4_FOR_PERIMETER
                            nr, nc = pix_r + dr, pix_c + dc
                            if !((1 <= nr <= rows) && (1 <= nc <= cols) && ((nr, nc) in pixel_set) )
                                obj_perimeter_edges += 1
                                is_on_boundary_for_plot = true # Pixel is on boundary if any 4-neighbor is background
                            end
                        end
                        if is_on_boundary_for_plot
                             # For plotting, store as (x,y)
                            push!(boundary_pixel_coords, Point2f(pix_c, pix_r))
                        end
                    end

                    obj_perimeter = Float32(obj_perimeter_edges)
                    obj_circularity = (obj_perimeter == 0) ? 0.0f0 : Float32((4 * π * obj_area) / (obj_perimeter^2))

                    push!(object_features_list, WormFeatures(
                        Point2f(centroid_x, centroid_y),
                        obj_area,
                        obj_perimeter,
                        obj_circularity,
                        current_object_pixels_tuples, # Store all pixels for potential detailed outline
                        boundary_pixel_coords # Store boundary pixels for simpler outline plotting
                    ))
                end
            end
        end
    end
    return object_features_list
end

"""
Processes a single image frame to find features of potential worms.
Returns a vector of WormFeatures.
"""
function find_worm_features_in_frame(img::AbstractMatrix{<:Gray}, min_area::Int, threshold::Float64)
    binary_mask = img .> threshold
    return find_objects_features(binary_mask, min_area)
end

"""
Main function to load TIFF stack, track multiple worms, and display paths & features.
"""
function track_worms_and_features_from_tiff(filepath::String)
    println("Loading TIFF stack from: $filepath")
    img_stack_data = TiffImages.load(filepath)

    frames_iterable = AbstractMatrix{<:Gray}[]
    num_frames = 0
    img_rows, img_cols = 0, 0

    if ndims(img_stack_data) == 3 && size(img_stack_data, 3) > 1 # Typical TIFF stack
        num_frames = size(img_stack_data, 3)
        img_rows, img_cols = size(img_stack_data, 1), size(img_stack_data, 2)
        frames_iterable = [Gray.(img_stack_data[:, :, i]) for i in 1:num_frames]
    elseif isa(img_stack_data, Vector) && !isempty(img_stack_data) # Vector of images
        matrix_elements = filter(x -> isa(x, AbstractMatrix), img_stack_data)
        num_frames = length(matrix_elements)
        if num_frames > 0
            img_rows, img_cols = size(matrix_elements[1],1), size(matrix_elements[1],2)
            frames_iterable = [Gray.(m) for m in matrix_elements]
        end
    else # Single image
        println("Warning: TIFF does not appear to be a stack or a list of frames. Processing as a single frame.")
        num_frames = 1
        img_rows, img_cols = size(img_stack_data,1), size(img_stack_data,2)
        frames_iterable = [Gray.(img_stack_data)]
    end

    if num_frames == 0
        error("No frames found in the TIFF file.")
    end
    println("Found $num_frames frames with dimensions $img_rows x $img_cols.")

    # --- Unit Conversion Factors ---
    cm_per_pixel = 1.0f0 / PIXELS_PER_CM
    area_cm2_per_pixel2 = cm_per_pixel^2
    seconds_per_frame = MSEC_PER_FRAME / 1000.0f0

    # --- Tracking Data Structures ---
    # Store last known features for active tracks
    active_tracks = Dict{Int, WormFeatures}()
    active_tracks_last_pos3d = Dict{Int, Point3f}() # For linking

    # Store completed track data
    completed_positions = Dict{Int, Vector{Point3f}}()
    completed_areas_px = Dict{Int, Vector{Int}}()
    completed_circularities = Dict{Int, Vector{Float32}}()
    completed_outlines_px = Dict{Int, Vector{Vector{Point2f}}}() # Store boundary pixels for each frame

    next_worm_id = 1

    println("Starting worm tracking...")
    for frame_idx in 1:num_frames
        current_frame_img = frames_iterable[frame_idx]
        current_features_list = find_worm_features_in_frame(current_frame_img, MIN_WORM_AREA, THRESHOLD_VALUE)

        num_current_centroids = length(current_features_list)
        matched_to_active_track = falses(num_current_centroids)
        new_active_tracks_this_frame = Dict{Int, WormFeatures}()
        new_active_tracks_last_pos3d_this_frame = Dict{Int, Point3f}()

        # Link active tracks
        for (id, last_pos3d) in active_tracks_last_pos3d
            min_dist = MAX_LINKING_DISTANCE
            best_match_idx = -1
            last_pos_2d = Point2f(last_pos3d[1], last_pos3d[2]) # Centroid X, Y in pixels

            for k in 1:num_current_centroids
                if !matched_to_active_track[k]
                    current_centroid_2d = current_features_list[k].centroid
                    dist = distance_2d(last_pos_2d, current_centroid_2d)
                    if dist < min_dist
                        min_dist = dist
                        best_match_idx = k
                    end
                end
            end

            if best_match_idx != -1
                matched_feature = current_features_list[best_match_idx]
                matched_point_3d = Point3f(matched_feature.centroid[1], matched_feature.centroid[2], frame_idx)

                push!(completed_positions[id], matched_point_3d)
                push!(completed_areas_px[id], matched_feature.area)
                push!(completed_circularities[id], matched_feature.circularity)
                push!(completed_outlines_px[id], matched_feature.boundary_pixels)


                new_active_tracks_this_frame[id] = matched_feature
                new_active_tracks_last_pos3d_this_frame[id] = matched_point_3d
                matched_to_active_track[best_match_idx] = true
            end
        end

        # Add new tracks for unmatched centroids
        for k in 1:num_current_centroids
            if !matched_to_active_track[k]
                new_id = next_worm_id
                next_worm_id += 1
                new_feature = current_features_list[k]
                new_point_3d = Point3f(new_feature.centroid[1], new_feature.centroid[2], frame_idx)

                completed_positions[new_id] = [new_point_3d]
                completed_areas_px[new_id] = [new_feature.area]
                completed_circularities[new_id] = [new_feature.circularity]
                completed_outlines_px[new_id] = [new_feature.boundary_pixels]


                new_active_tracks_this_frame[new_id] = new_feature
                new_active_tracks_last_pos3d_this_frame[new_id] = new_point_3d
            end
        end
        active_tracks = new_active_tracks_this_frame
        active_tracks_last_pos3d = new_active_tracks_last_pos3d_this_frame

        if frame_idx % 50 == 0 || frame_idx == num_frames
             println("Processed frame $frame_idx / $num_frames")
        end
    end

    println("Tracking complete. Found $(length(completed_positions)) track segments.")

    # --- Plotting ---
    GLMakie.set_theme!(theme_light()) # Or theme_dark()
    fig = Figure(size=(1600, 1200))
    ga = fig[1, 1] = GridLayout()
    gb = fig[2, 1:2] = GridLayout() # Span two columns for feature plots

    ax_tracks = Axis3(ga[1,1],
                      title="Worm Trajectories",
                      xlabel="X (cm)", ylabel="Y (cm)", zlabel="Time (s)",
                      aspect=(1, img_rows/img_cols * cm_per_pixel / cm_per_pixel, 0.5) # Adjust aspect for better visualization
                     )

    ax_outline_example = Axis(ga[1,2], title="Example Worm Outline", aspect=DataAspect())
    ax_size = Axis(gb[1,1], xlabel="Time (s)", ylabel="Area (cm²)", title="Worm Size Over Time")
    ax_circularity = Axis(gb[1,2], xlabel="Time (s)", ylabel="Circularity", title="Worm Circularity Over Time")

    # Get IDs of tracks long enough to be plotted
    min_track_length = 10 # Minimum frames for a track to be considered for feature plotting
    track_ids_to_plot = [id for (id, track) in completed_positions if length(track) >= min_track_length]

    if isempty(track_ids_to_plot)
        println("No tracks long enough to plot features.")
    else
        # Generate distinguishable colors for tracks
        plot_colors = distinguishable_colors(length(track_ids_to_plot), [RGB(1,1,1), RGB(0,0,0)], lchoices=range(20, stop=70, length=15))

        for (plot_idx, id) in enumerate(track_ids_to_plot)
            track_points_px = completed_positions[id]
            track_areas_px = completed_areas_px[id]
            track_circularities = completed_circularities[id]

            # Convert to physical units
            track_points_cm_s = [Point3f(p[1]*cm_per_pixel, (img_rows - p[2])*cm_per_pixel, p[3]*seconds_per_frame) for p in track_points_px]
            track_times_s = [p[3]*seconds_per_frame for p in track_points_px]
            track_areas_cm2 = [area * area_cm2_per_pixel2 for area in track_areas_px]

            color_to_use = plot_colors[plot_idx]

            lines!(ax_tracks, track_points_cm_s, color=color_to_use, linewidth=2, label="Worm $id")

            lines!(ax_size, track_times_s, track_areas_cm2, color=color_to_use, linewidth=1.5, label="Worm $id")
            scatter!(ax_size, track_times_s, track_areas_cm2, color=color_to_use, markersize=4)

            lines!(ax_circularity, track_times_s, track_circularities, color=color_to_use, linewidth=1.5, label="Worm $id")
            scatter!(ax_circularity, track_times_s, track_circularities, color=color_to_use, markersize=4)
        end
        # Add legends to feature plots
        Legend(gb[1,3], ax_size, "Worm ID", framevisible=false) # A single legend for both
    end


    # Plot outline example for the first "good" track in the middle of the video
    if !isempty(track_ids_to_plot)
        example_track_id = track_ids_to_plot[1]
        track_len = length(completed_positions[example_track_id])
        middle_frame_in_track_idx = div(track_len, 2) + 1

        if middle_frame_in_track_idx > 0 && middle_frame_in_track_idx <= length(completed_outlines_px[example_track_id])
            original_frame_idx = Int(completed_positions[example_track_id][middle_frame_in_track_idx][3]) # Get original frame index

            # Display the original binary image for context
            binary_middle_frame = Gray.(frames_iterable[original_frame_idx]) .> THRESHOLD_VALUE
            # Convert image to cm scale for axes
            x_coords_cm = (1:img_cols) .* cm_per_pixel
            y_coords_cm = (1:img_rows) .* cm_per_pixel
            image!(ax_outline_example, x_coords_cm, y_coords_cm, rotr90(binary_middle_frame), colormap=:grays)

            # Plot outline
            outline_pixels_px = completed_outlines_px[example_track_id][middle_frame_in_track_idx]
            outline_points_cm = [Point2f(p[1]*cm_per_pixel, (img_rows - p[2])*cm_per_pixel) for p in outline_pixels_px] # Correct Y for image origin
            if !isempty(outline_points_cm)
                scatter!(ax_outline_example, outline_points_cm, color=plot_colors[1], markersize=3, label="Outline")
            end

            # Plot centroid
            centroid_px = Point2f(completed_positions[example_track_id][middle_frame_in_track_idx][1], completed_positions[example_track_id][middle_frame_in_track_idx][2])
            centroid_cm = Point2f(centroid_px[1]*cm_per_pixel, (img_rows - centroid_px[2])*cm_per_pixel)
            scatter!(ax_outline_example, centroid_cm, color=:red, marker=:cross, markersize=15, label="Centroid")
            ax_outline_example.title = "Outline (Worm $example_track_id, Frame $original_frame_idx, Time $(round(original_frame_idx*seconds_per_frame, digits=2))s)"
        end
    else
        text!(ax_outline_example, "No suitable tracks for outline example.", position = (0.1, 0.5), align = (:left, :center))
    end
    xlims!(ax_outline_example, 0, img_cols * cm_per_pixel)
    ylims!(ax_outline_example, 0, img_rows * cm_per_pixel)


    display(fig)
    println("Displayed plots. Ensure GLMakie window is interactive.")
end


# --- Functions for Video Comparison (TODO: Implement fully in a separate context or script) ---
function _internal_track_worms_from_tiff_for_diff(filepath::String)
    println("Placeholder: Running simplified tracking for difference calculation on $filepath")
    # This would be a simplified version of track_worms_and_features_from_tiff
    # It should return a Dict{Int, Vector{Point3f}} of pixel coordinates and frame indices
    # For now, it's just a placeholder.
    # To make this work, you'd essentially call the core tracking logic from track_worms_and_features_from_tiff
    # and extract `completed_positions`.
    # Example (conceptual, needs actual implementation of the tracking part):
    # global MIN_WORM_AREA, THRESHOLD_VALUE, MAX_LINKING_DISTANCE (if not const or passed as args)
    # ... (load image stack as in the main function) ...
    # ... (tracking loop as in the main function, populating a local completed_positions) ...
    # return local_completed_positions
    return Dict{Int,Vector{Point3f}}() # Return empty for placeholder
end

function plot_path_difference(filepath1::String, filepath2::String)
    println("Attempting to plot path difference between $filepath1 and $filepath2")
    # For this to work, _internal_track_worms_from_tiff_for_diff needs to be implemented
    # to perform tracking and return the completed_positions.
    raw_path1_map = _internal_track_worms_from_tiff_for_diff(filepath1)
    raw_path2_map = _internal_track_worms_from_tiff_for_diff(filepath2)

    # --- Unit Conversion Factors (repeated for standalone use if necessary) ---
    cm_per_pixel_val = 1.0f0 / PIXELS_PER_CM
    seconds_per_frame_val = MSEC_PER_FRAME / 1000.0f0

    # Assume one primary worm track, or merge/select one (e.g., longest)
    path1_px = isempty(raw_path1_map) ? Point3f[] : raw_path1_map[argmax(length, raw_path1_map)]
    path2_px = isempty(raw_path2_map) ? Point3f[] : raw_path2_map[argmax(length, raw_path2_map)]

    if isempty(path1_px) || isempty(path2_px)
        println("One or both videos resulted in no tracks for difference plot.")
        return
    end

    println("Path 1 length (frames): $(length(path1_px)), Path 2 length (frames): $(length(path2_px))")

    # Convert to cm and seconds
    path1 = [Point3f(p[1]*cm_per_pixel_val, p[2]*cm_per_pixel_val, p[3]*seconds_per_frame_val) for p in path1_px]
    path2 = [Point3f(p[1]*cm_per_pixel_val, p[2]*cm_per_pixel_val, p[3]*seconds_per_frame_val) for p in path2_px]

    min_len = min(length(path1), length(path2))
    if min_len == 0
        println("Tracks are empty after processing. Cannot compute difference.")
        return
    end

    differences_cm = Float32[]
    times_s = Float32[]

    # Assuming frame indices in Point3f's z component correspond and can be used for time alignment.
    # For a more robust solution, time-based interpolation might be needed if frame rates differ or tracks are sparse.
    map1_time_to_pos = Dict(p[3] => Point2f(p[1], p[2]) for p in path1)
    map2_time_to_pos = Dict(p[3] => Point2f(p[1], p[2]) for p in path2)

    common_times = sort(collect(intersect(keys(map1_time_to_pos), keys(map2_time_to_pos))))

    if isempty(common_times)
        println("No common time points found between tracks based on frame indices. Cannot compute difference without interpolation.")
        return
    end
    println("Found $(length(common_times)) common time points for difference calculation.")

    for t in common_times
        pos1 = map1_time_to_pos[t]
        pos2 = map2_time_to_pos[t]
        diff_vec = pos1 - pos2 # Difference in (x,y) plane in cm
        push!(differences_cm, norm(diff_vec))
        push!(times_s, t) # Time is already in seconds
    end

    if isempty(differences_cm)
        println("No differences calculated. Check track alignment and data.")
        return
    end

    fig_diff = Figure(size=(800,600))
    ax_diff = Axis(fig_diff[1,1], xlabel="Time (s)", ylabel="Path Difference (cm)", title="Difference between Worm Paths")
    lines!(ax_diff, times_s, differences_cm, color=:purple, linewidth=2)
    scatter!(ax_diff, times_s, differences_cm, color=:purple, markersize=4)
    display(fig_diff)
    println("Displayed path difference plot.")
end


# --- Main Execution ---
track_worms_and_features_from_tiff(TIFF_FILE_PATH)

# --- TODO: Video Difference Plot ---
# This part requires two separate TIFF files and a fully implemented
# _internal_track_worms_from_tiff_for_diff or a similar function.
# For now, this is a placeholder and would need to be called with actual file paths.
# println("\n--- Video Difference Plot (TODO) ---")
# const TIFF_FILE_PATH_2 = "/path/to/your/second_video.tif" # !!! REPLACE THIS !!!
# if TIFF_FILE_PATH != TIFF_FILE_PATH_2 && isfile(TIFF_FILE_PATH_2)
#     plot_path_difference(TIFF_FILE_PATH, TIFF_FILE_PATH_2)
# else
#     println("To run path difference plot, define TIFF_FILE_PATH_2 with a different valid video.")
# end