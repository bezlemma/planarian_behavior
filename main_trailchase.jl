# Main script for processing and plotting worm tracking
# Usage: define filepaths, thresholds, and plot options, then call `track_and_plot`.

include("Stab_Interp_Track_Worms.jl")
using .Stab_Interp_Track_Worms: track_worms

# Define input files and parameters
filepaths = [
    "/Users/dl0346/Documents/PlanarianVideos/May28/TrailFollow/ATM4_1.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/TrailFollow/ATM4_2.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/TrailFollow/ATM4_3.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/TrailFollow/ATM4_4.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/TrailFollow/ATM4_5.tif",
]

thresholds = [1.5, 1.5, 1.5, 1.5, 1.5]  # Per-file thresholds

# Call the processing function
results, stacks = track_worms(filepaths; thresholds=thresholds)

# Plot 3D trajectories
include("Plotting_Aux.jl")
using .Plotting_Aux: plot_3d_trajectory
plot_3d_trajectory(results)

