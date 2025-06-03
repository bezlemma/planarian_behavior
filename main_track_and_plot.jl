# Main script for processing and plotting worm tracking
# Usage: define filepaths, thresholds, and plot options, then call `track_and_plot`.

include("Stab_Interp_Track_Worms_LightDark.jl")
using .Stab_Interp_Track_Worms_LightDark: track_worms

# Define input files and parameters
filepaths = [
    "/Users/dl0346/Documents/PlanarianVideos/May28/HalfRed/ATM1_1.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/HalfRed/ATM1_2.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/HalfRed/ATM1_3.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/HalfRed/ATM1_4.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/HalfRed/ATM1_5.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/HalfRed/ATM1_6.tif",
]

thresholds = [1.1, 1.1, 1.1, 1.1, 1.1, 1.1]  # Per-file thresholds

# Call the processing function
results, stacks = track_worms(filepaths; thresholds=thresholds)

# Generate light/dark plot if requested
include("Plotting_Aux.jl")
using .Plotting_Aux: plot_lightdark


plot_lightdark(results)

