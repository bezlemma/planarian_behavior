# Main script for processing and plotting worm tracking
# Usage: define filepaths, thresholds, and plot options, then call `track_and_plot`.

include("Stab_Interp_Track_Worms.jl")
using .Stab_Interp_Track_Worms: track_worms

# Define input files and parameters
filepaths = [
    "/Users/dl0346/Documents/PlanarianVideos/May28/HalfRed/ATM1_1.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/HalfRed/ATM1_2.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/HalfRed/ATM1_3.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/HalfRed/ATM1_4.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/HalfRed/ATM1_5.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/HalfRed/ATM1_6.tif",

    "/Users/dl0346/Documents/PlanarianVideos/May28/HalfRed/ATMUNC1_1.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/HalfRed/ATMUNC1_2.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/HalfRed/ATMUNC1_3.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/HalfRed/ATMUNC1_4.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/HalfRed/ATMUNC1_5.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/HalfRed/ATMUNC1_6.tif",

    "/Users/dl0346/Documents/PlanarianVideos/May28/HalfRed/ATM4_1.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/HalfRed/ATM4_2.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/HalfRed/ATM4_3.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/HalfRed/ATM4_4.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/HalfRed/ATM4_5.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/HalfRed/ATM4_6.tif",

    "/Users/dl0346/Documents/PlanarianVideos/May28/HalfRed/ATMUNC4_1.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/HalfRed/ATMUNC4_2.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/HalfRed/ATMUNC4_3.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/HalfRed/ATMUNC4_4.tif", #TODO: This isn't tracked well enough
    "/Users/dl0346/Documents/PlanarianVideos/May28/HalfRed/ATMUNC4_5.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/HalfRed/ATMUNC4_6.tif",

    "/Users/dl0346/Documents/PlanarianVideos/May28/HalfRed/ATM6_1.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/HalfRed/ATM6_2.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/HalfRed/ATM6_3.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/HalfRed/ATM6_4.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/HalfRed/ATM6_5.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/HalfRed/ATM6_6.tif",

    "/Users/dl0346/Documents/PlanarianVideos/May28/HalfRed/NDK1_1.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/HalfRed/NDK1_2.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/HalfRed/NDK1_3.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/HalfRed/NDK1_4.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/HalfRed/NDK1_5.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/HalfRed/NDK1_6.tif",

    "/Users/dl0346/Documents/PlanarianVideos/May28/HalfRed/NDKUNC1_1.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/HalfRed/NDKUNC1_2.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/HalfRed/NDKUNC1_3.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/HalfRed/NDKUNC1_4.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/HalfRed/NDKUNC1_5.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/HalfRed/NDKUNC1_6.tif",

    "/Users/dl0346/Documents/PlanarianVideos/May28/HalfRed/TEC1_1.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/HalfRed/TEC1_2.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/HalfRed/TEC1_3.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/HalfRed/TEC1_4.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/HalfRed/TEC1_5.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/HalfRed/TEC1_6.tif",

    "/Users/dl0346/Documents/PlanarianVideos/May28/HalfRed/TECUNC1_1.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/HalfRed/TECUNC1_2.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/HalfRed/TECUNC1_3.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/HalfRed/TECUNC1_4.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/HalfRed/TECUNC1_5.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/HalfRed/TECUNC1_6.tif",
]

thresholds = [1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 
                1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 
                1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 
                1.1, 1.1, 1.1, 1.1, 1.1, 1.1,

                1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 
                1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 
                1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 
                1.1, 1.1, 1.1, 1.1, 1.1, 1.1]  # Per-file thresholds

# Call the processing function
results, stacks = track_worms(filepaths; thresholds=thresholds)


## TODO: It's not clear what to do here, as in the current data I don't see worms eating the food.