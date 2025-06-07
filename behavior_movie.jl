# Main script for processing and plotting worm tracking
# Usage: define filepaths, thresholds, and plot options, then call `track_and_plot`.

include("Stab_Interp_Track_Worms.jl")
using .Stab_Interp_Track_Worms: track_worms

# Define input files and parameters
filepath= ["/Users/dl0346/Documents/PlanarianVideos/May28/TrailFollow/ATM1_3.tif"]
threshold = [0.05] # Represents 0.05% of brightest pixels to be selected


# Call the processing function
results, stacks = track_worms(filepath; pixel_percentages=threshold)

# Plots
include("Plotting_Aux.jl")
using .Plotting_Aux

#
using GLMakie
using ColorSchemes
using Colors

# Automated movie recording of worm behavior
using .Plotting_Aux: behavior_colors
# Prepare data
stack = stacks[end]
worm = results[end]
# Map frame index to centroid, area, behavior
worm_info = Dict{Int, Tuple{Point2f, Int, Symbol}}()
for i in 1:length(worm.positions_)
    p3 = worm.positions_[i]
    f_idx = Int(round(p3[3]))
    cent = Point2f(p3[2], p3[1])
    worm_info[f_idx] = (cent, worm.areas_2[i], worm.behaviors[i])
end

# Build figure and plots
fig = Figure(size = (600,600))
ax = Axis(fig[1,1], aspect = DataAspect())
# Initialize image observable for dynamic frame updates
img_obs = Observable(Gray.(stack[:,:,1]))
image!(ax, img_obs, colormap = :grays, colorrange = (0,1))
pos_obs = Observable(Point2f(NaN, NaN))
col_obs = Observable(RGBA(1,1,1,0.7))
size_obs = Observable(50.0)
scatter!(ax,
    lift(pos_obs) do p; p[1] end,
    lift(pos_obs) do p; p[2] end,
    color = col_obs,
    markersize = size_obs,
    strokecolor = :transparent)

# Record frames into MP4
nframes = size(stack,3)
record(fig, "behavior_movie.mp4", 1:nframes; framerate = 100) do i
    # update image for current frame
    img_obs[] = Gray.(stack[:,:,i])
    if haskey(worm_info, i)
        cent, area, beh = worm_info[i]
        pos_obs[] = cent
        col_obs[] = behavior_colors[beh]
        size_obs[] = clamp(sqrt(max(0.0, Float64(area)) / π), 2, 25) * 10
    else
        pos_obs[] = Point2f(NaN, NaN)
    end
end
