# Main script for processing and plotting worm tracking
# Usage: define filepaths, thresholds, and plot options, then call `track_and_plot`.

include("Stab_Interp_Track_Worms.jl")
using .Stab_Interp_Track_Worms: track_worms

# Define input files and parameters
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

thresholds = [1.5, 1.5, 1.5, 1.5, 1.5,
1.5, 1.5, 1.5, 1.5, 1.5,
1.5, 1.5, 1.5, 1.5, 1.5,
1.5, 1.5, 1.5, 1.5, 1.5,

1.5, 1.5, 1.5, 1.5, 1.5,
1.5, 1.5, 1.5, 1.5, 1.5,
1.5, 1.5, 1.5, 1.5, 1.5,
1.5, 1.5, 1.5, 1.5, 1.5,
]  # Per-file thresholds

# Call the processing function
results, stacks = track_worms(filepaths; thresholds=thresholds)

# Plots
include("Plotting_Aux.jl")
using .Plotting_Aux

##
using GLMakie
using ColorSchemes

cond_labels = ["ATM1", "ATMUNC1", "ATM4", "ATMUNC4", "TEC1", "TECUNC1", "NDK1", "NDKUNC1"]
plot_behavior_distribution(results, cond_labels)


##
# plot cumulative trail chase distances with color gradients
function plot_trail_chase_cumulative(results)
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

    p1 = lines!(ax, ts, ATM1; color=ts, colormap=:blues, linewidth=5)
    p2 = lines!(ax, ts, ATMUNC1; color=:blue, linewidth=5)
    p3 = lines!(ax, ts, ATM4; color=ts, colormap=:greens, linewidth=5)
    p4 = lines!(ax, ts, ATMUNC4; color=:green, linewidth=5)
    p5 = lines!(ax, ts, NDK1; color=ts, colormap=:reds, linewidth=5)
    p6 = lines!(ax, ts, NDKUNC1; color=:red, linewidth=5)
    p7 = lines!(ax, ts, TEC1; color=ts, colormap=ColorSchemes.Oranges, linewidth=5)
    p8 = lines!(ax, ts, TECUNC1; colormap=:orange, linewidth=5)

    mid_blue = Makie.to_colormap(:blues)[div(length(Makie.to_colormap(:blues)), 2)]
    mid_red = Makie.to_colormap(:reds)[div(length(Makie.to_colormap(:reds)), 2)]
    mid_green = Makie.to_colormap(:greens)[div(length(Makie.to_colormap(:greens)), 2)]
    mid_orange = ColorSchemes.Oranges[0.5]

    legend_elements = [
        LineElement(color = mid_blue, linewidth = 5),
        LineElement(color = :blue, linewidth = 5),
        LineElement(color = mid_green, linewidth = 5),
        LineElement(color = :green, linewidth = 5),
        LineElement(color = mid_red, linewidth = 5),
        LineElement(color = :red, linewidth = 5),
        LineElement(color = mid_orange, linewidth = 5),
        LineElement(color = :orange, linewidth = 5)
    ]
    
    # Create the legend using the custom elements
    axislegend(ax, legend_elements, ["ATM1", "ATMUNC1", "ATM4", "ATMUNC4","NDK1", "NDKUNC1", "TEC1", "TECUNC1"]; position = :lt)


    display(fig)
    return fig
end

# Generate the final plot
plot_trail_chase_cumulative(results)
