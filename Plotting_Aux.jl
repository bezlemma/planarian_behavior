module Plotting_Aux

using GLMakie, GeometryBasics, Colors

# Unit conversion constants
const CM_PER_PIXEL = 1.0 / 58.6
const MSEC_PER_FRAME = 50.0 * 2.0

export plot_all_distances_from_center, calculate_msd, plot_all_msd, plot_sum_pairwise_differences_over_time, view_stack_and_worm, plot_lightdark

function plot_all_distances_from_center(all_worm_data)
    fig = Figure(size=(1000, 750))
    ax = Axis(fig[1,1], xlabel="Time (s)", ylabel="Distance from Center (cm)", title="Distance from Image Center Over Time")
    colors = distinguishable_colors(max(1,length(all_worm_data)), [RGB(1,1,1), RGB(0,0,0)], lchoices=range(20, stop=70, length=15))
    for (i, data) in enumerate(all_worm_data)
        center_x_px = data.img_cols / 2.0; center_y_px = data.img_rows / 2.0
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

function calculate_msd(positions, times)
    n_points = length(positions)
    if n_points < 2 || length(times) != n_points return Dict{Float32, Float32}() end
    msd_values = Dict{Float32, Vector{Float32}}()
    for lag_frames in 1:n_points-1
        for i in 1:(n_points - lag_frames)
            p1 = positions[i]; p2 = positions[i + lag_frames]
            diff_vec = p2 - p1; dist_sq = sum(diff_vec.^2)
            time_lag_s = abs(times[i + lag_frames] - times[i])
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

function plot_all_msd(all_worm_data)
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
        actual_min_lag = max(min_lag_overall, 1e-3); actual_max_lag = max_lag_overall
        if actual_min_lag < actual_max_lag
            ref_lags = 10 .^ LinRange(log10(actual_min_lag), log10(actual_max_lag), 50)
            anchor_y_at_min_lag = isfinite(min_msd_val_overall) ? max(min_msd_val_overall, 1e-9) : 1e-4
            y1_vals = anchor_y_at_min_lag .* (ref_lags ./ ref_lags[1]).^1.0; lines!(ax, ref_lags, y1_vals, color=RGB(0.5,0.5,0.5), linestyle=:dash, label="α=1")
            y2_vals = anchor_y_at_min_lag .* (ref_lags ./ ref_lags[1]).^2.0; lines!(ax, ref_lags, y2_vals, color=RGB(0.7,0.7,0.7), linestyle=:dot, label="α=2")
            y05_vals = anchor_y_at_min_lag .* (ref_lags ./ ref_lags[1]).^0.5; lines!(ax, ref_lags, y05_vals, color=RGB(0.3,0.3,0.3), linestyle=:dashdot, label="α=0.5")
        end
    end
    if plotted_anything axislegend(ax, position=:rb) end
    display(GLMakie.Screen(), fig); return fig
end

function plot_sum_pairwise_differences_over_time(all_worm_data)
    fig = Figure(size=(1000, 750))
    ax = Axis(fig[1,1], xlabel="Time (s)", ylabel="Sum of Pairwise Distances (cm)", title="Sum of Pairwise Worm Distances Over Time")
    if length(all_worm_data) < 2 display(GLMakie.Screen(), fig); return fig end
    all_frame_indices_with_data = Set{Int}(); frame_to_time_map = Dict{Int, Float32}()
    worm_pos_at_frame_idx = [Dict{Int, Point2f}() for _ in 1:length(all_worm_data)]
    for (data_idx, data) in enumerate(all_worm_data)
        for p3d_px in data.positions_px
            frame_idx = round(Int, p3d_px[3]); time_s = frame_idx * MSEC_PER_FRAME / 1000.0
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
            current_sum_pairwise_dist = 0.0
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



function view_stack_and_worm(binary_stack_to_show,tracked_worm_data)
    marker_color=RGBAf(1.0, 0.0, 0.0, 0.7) 
    rows, cols, num_frames = size(binary_stack_to_show)
    worm_info_for_frame = Dict{Int, Tuple{Point2f, Int}}() 
    for i in 1:length(tracked_worm_data.positions_)
        pt3 = tracked_worm_data.positions_[i]
        frame_idx = Int(round(pt3[3]))
        centroid = Point2f(pt3[2], pt3[1]) # (col_, row_) for GLMakie point
        area = tracked_worm_data.areas_2[i]
        worm_info_for_frame[frame_idx] = (centroid, area)
    end
    fig = Figure(size = (cols > rows ? (800, 800 * rows/cols + 100) : (800 * cols/rows + 100, 800)))
    ax_img = Axis(fig[1, 1])
    slider = Slider(fig[2, 1], range=1:num_frames, startvalue=1)
    frame_idx_obs = slider.value
    img_slice_obs = lift(frame_idx_obs) do f_idx
        return Gray.(view(binary_stack_to_show, :, :, f_idx)) 
    end
    image!(ax_img, img_slice_obs, interpolate=false, colormap=:grays, colorrange=(0,1))
    worm_marker_obs = Observable([Circle(Point2f(NaN, NaN), 0)])
    on(frame_idx_obs) do f_idx
        if haskey(worm_info_for_frame, f_idx)
            centroid, wormarea = worm_info_for_frame[f_idx]
            radius = sqrt(max(0.0, Float64(wormarea)) / π)
            display_radius = clamp(Float32(radius), 2, 25)
            worm_marker_obs[] = [Circle(centroid, display_radius)]
        else
            worm_marker_obs[] = [Circle(Point2f(NaN, NaN), 0)]
        end
    end
    poly!(ax_img, worm_marker_obs, color=marker_color, strokecolor=:transparent)
    display(GLMakie.Screen(), fig)
    return fig
end

# Function to plot light/dark trajectories
function plot_lightdark(worm_results)
    fig = Figure(size=(1000, 1000))
    img_width_cm = 7.5; img_height_cm = 7.5

    ax = Axis(fig[1,1], xlabel="X Position (cm)", ylabel="Y Position (cm)",
        xlabelsize=14, ylabelsize=14, xgridvisible=false,
        ygridvisible=false, limits=((0, img_width_cm), (0, img_height_cm)))

    # Draw boundary circle
    centerX = img_width_cm / 2.0; centerY = img_height_cm / 2.0
    circle_radius = img_height_cm / 2.05
    circle_geom = Circle(Point2f(centerX, centerY), circle_radius)
    poly!(ax, circle_geom, color=:transparent, strokecolor=:black, strokewidth=2)

    # Shade bottom half
    num_points_arc = 100
    theta_arc = range(pi, 0, length=num_points_arc)
    arc_x = centerX .+ circle_radius .* cos.(theta_arc)
    arc_y_bottom = centerY .- circle_radius .* sin.(theta_arc)
    half_circle_vertices_bottom = [Point2f(x, y) for (x, y) in zip(arc_x, arc_y_bottom)]
    poly!(ax, half_circle_vertices_bottom, color=RGBA(0.6, 0.1, 0.1, 0.8))

    # Plot trajectories
    for worm in worm_results
        points = worm.positions_cm
        x_pos = [pt[1] for pt in points]
        y_pos = [pt[2] for pt in points]
        subsample_ix = 1:2:length(points)
        x_sub = x_pos[subsample_ix]; y_sub = y_pos[subsample_ix]
        custom_cmap = cgrad([RGBA(colorant"cyan", 0.1), RGBA(colorant"blue", 0.1)])
        if !isempty(x_sub)
            scatter!(ax, x_sub, y_sub, color=y_sub, colormap=custom_cmap, markersize=15)
        end
    end

    display(GLMakie.Screen(), fig)
    return fig
end

end # module Plotting_Aux