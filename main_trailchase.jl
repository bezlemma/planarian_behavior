# Main script for processing and plotting worm tracking
# Usage: define filepaths, thresholds, and plot options, then call `track_and_plot`.

include("Stab_Interp_Track_Worms.jl")
using .Stab_Interp_Track_Worms: track_worms

# Define input files and parameters
filepaths_test = ["/Users/dl0346/Documents/PlanarianVideos/May28/TrailFollow/ATMUNC4_3.tif"]
thresholds_test = [0.05] # Represents 0.05% of brightest pixels to be selected

filepaths = [
    "/Users/dl0346/Documents/PlanarianVideos/May28/TrailFollow/ATM1_1.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/TrailFollow/ATM1_2.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/TrailFollow/ATM1_3.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/TrailFollow/ATM1_4.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/TrailFollow/ATM1_5.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/TrailFollow/ATMUNC1_1.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/TrailFollow/ATMUNC1_2.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/TrailFollow/ATMUNC1_3.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/TrailFollow/ATMUNC1_4.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/TrailFollow/ATMUNC1_5.tif",

    "/Users/dl0346/Documents/PlanarianVideos/May28/TrailFollow/ATM4_1.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/TrailFollow/ATM4_2.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/TrailFollow/ATM4_3.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/TrailFollow/ATM4_4.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/TrailFollow/ATM4_5.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/TrailFollow/ATMUNC4_1.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/TrailFollow/ATMUNC4_2.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/TrailFollow/ATMUNC4_3.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/TrailFollow/ATMUNC4_4.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/TrailFollow/ATMUNC4_5.tif",

     "/Users/dl0346/Documents/PlanarianVideos/May28/TrailFollow/ATM6_1.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/TrailFollow/ATM6_2.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/TrailFollow/ATM6_3.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/TrailFollow/ATM6_4.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/TrailFollow/ATM6_5.tif",

    "/Users/dl0346/Documents/PlanarianVideos/May28/TrailFollow/NDK1_1.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/TrailFollow/NDK1_2.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/TrailFollow/NDK1_3.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/TrailFollow/NDK1_4.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/TrailFollow/NDK1_5.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/TrailFollow/NDKUNC1_1.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/TrailFollow/NDKUNC1_2.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/TrailFollow/NDKUNC1_3.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/TrailFollow/NDKUNC1_4.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/TrailFollow/NDKUNC1_5.tif",

    "/Users/dl0346/Documents/PlanarianVideos/May28/TrailFollow/TEC1_1.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/TrailFollow/TEC1_2.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/TrailFollow/TEC1_3.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/TrailFollow/TEC1_4.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/TrailFollow/TEC1_5.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/TrailFollow/TECUNC1_1.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/TrailFollow/TECUNC1_2.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/TrailFollow/TECUNC1_3.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/TrailFollow/TECUNC1_4.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/TrailFollow/TECUNC1_5.tif",
]

thresholds = [0.05, 0.05, 0.05, 0.05, 0.05,   #ATM1
              0.05, 0.05, 0.05, 0.05, 0.05, 
              0.05, 0.05, 0.05, 0.05, 0.05,   #ATM4
              0.05, 0.05, 0.05, 0.05, 0.05,   #ATMUNC4
            0.05, 0.05, 0.05, 0.05, 0.05,     #ATM6
            0.05, 0.05, 0.05, 0.05, 0.05,  
            0.05, 0.05, 0.05, 0.05, 0.05,    #NDKUNC
            0.05, 0.05, 0.05, 0.05, 0.05,  
            0.05, 0.05, 0.05, 0.05, 0.05] 

# Call the processing function
results, stacks = track_worms(filepaths; pixel_percentages=thresholds)

# Plots
include("Plotting_Aux.jl")
using .Plotting_Aux

#
using GLMakie
using ColorSchemes

cond_labels = [  "ATM1", "ATM4", "ATM6", "ATMUNC1", "ATMUNC4"] #Put whatever conditions you would like to plot here
plot_behavior_distribution(results, cond_labels)

view_stack_and_worm(stacks[end],results[end]) #Use the live-viewer

##
# plot cumulative trail chase distances with color gradients
function plot_trail_chase_diffs(results)
    ts = results[1].times_s
    # helper for cumulative series of min distances
    function cum_series(indices)
        positions = [results[i].positions_cm for i in indices]
        vals = zeros(Float32, length(ts))
        for k in 2:length(positions)
            len_k = length(positions[k])
            for t in 1:length(ts)
                if t > len_k continue end
                p = positions[k][t]
                min_d = minimum([
                    minimum([sqrt((p[1]-q[1])^2 + (p[2]-q[2])^2) for q in positions[m]])
                    for m in 1:(k-1)
                ])
                vals[t] += abs(min_d)
            end
        end
        return cumsum(vals)
    end
    # compute cumulative distances for each group
    ATM1 = cum_series(1:5)
    ATMUNC1 = cum_series(6:10)
    ATM4 = cum_series(11:15)
    ATMUNC4 = cum_series(16:20)
    NDK1 = cum_series(21:25)
    NDKUNC1 = cum_series(26:30)
    TEC1 = cum_series(31:35)
    TECUNC1 = cum_series(36:40)

    fig = Figure(size=(600, 400))
    ax = Axis(fig[1, 1], xlabel="Time (s)", ylabel="Accumulated Min Distance (cm)", title="Trail Chase Distances")

    lines!(ax, ts, abs.(ATM1 .- ATMUNC1); color=ts, colormap=ColorSchemes.Blues, linewidth=5)
    lines!(ax, ts, abs.(ATM4 .- ATMUNC4); color=ts, colormap=ColorSchemes.Greens, linewidth=5)
    lines!(ax, ts, abs.(NDK1 .- NDKUNC1); color=ts, colormap=ColorSchemes.Reds, linewidth=5)
    lines!(ax, ts, abs.(TEC1 .- TECUNC1); color=ts, colormap=ColorSchemes.Purples, linewidth=5)

    
    mid_blue = ColorSchemes.Blues[0.7]
    mid_red =  ColorSchemes.Reds[0.7]
    mid_green = ColorSchemes.Greens[0.7]
    mid_Purples = ColorSchemes.Purples[0.7]

    legend_elements = [
        LineElement(color = mid_blue, linewidth = 5),
        LineElement(color = mid_green, linewidth = 5),
        LineElement(color = mid_red, linewidth = 5),
        LineElement(color = mid_Purples, linewidth = 5),
    ]
    

    axislegend(ax, legend_elements, ["ATM1/ATMUNC1", "ATM4/ATMUNC4","NDK1/NDKUNC1", "TEC1/TECUNC1"]; position = :lt)

    display(fig)
    return fig
end

# Generate the final plot
plot_trail_chase_diffs(results)
