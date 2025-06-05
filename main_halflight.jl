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

# Generate light/dark plot if requested
#include("Plotting_Aux.jl")
#using .Plotting_Aux: plot_lightdark
#plot_lightdark(results)
##



# Compute and plot % time spent in dark (lower half) per condition
# Define conditions and grouping
conditions = ["ATM1", "ATMUNC1", "ATM4", "ATMUNC4", "NDK1", "NDKUNC1", "TEC1", "TECUNC1"]
n_reps = 6
n_conditions = length(conditions)
pct_dark = Float64[]
# collect raw replicate percentages
rep_pcts_list = Vector{Vector{Float64}}()

# Calculate mean % time in dark for each condition
enum_indices = 1:n_conditions
for (i, cond) in enumerate(conditions)
    idxs = ((i-1)*n_reps+1):(i*n_reps)
    rep_pcts = Float64[]
    for idx in idxs
        ys = [pt[2] for pt in results[idx].positions_cm]
        if isempty(ys)
            @warn "No position data for file $(results[idx].filepath), skipping."
            continue
        end
        thr = maximum(ys) / 2
        pct = count(y -> y < thr, ys) / length(ys) * 100
        push!(rep_pcts, pct)
    end
    # store replicate percentages
    push!(rep_pcts_list, rep_pcts)
    # compute and store mean
    if isempty(rep_pcts)
        @warn "No valid replicates for condition $cond, setting NaN."
        push!(pct_dark, NaN)
    else
        push!(pct_dark, mean(rep_pcts))
    end
end

# Plot bar graph with error bars and significance
using GLMakie, Statistics, HypothesisTests

function significance_level(p)
    if p < 0.001; println("p=$p, ***"); return "***"
    elseif p < 0.01; println("p=$p, **"); return "**"
    elseif p < 0.05; println("p=$p, *"); return "*"
    else println("p=$p, ns"); return "ns"; end
end

# compute standard error of the mean for each condition
sems = [ isempty(v) ? NaN : std(v) / sqrt(length(v)) for v in rep_pcts_list ]

fig_dark = Figure()
ax_dark = Axis(fig_dark[1,1], xticks=(enum_indices, conditions), ylabel="% Time in Dark",
    title="% Time in Dark by Condition", xticklabelrotation=45)

# set y limits to cover error bars and room for annotations
ylims!(ax_dark, 0, maximum(pct_dark .+ sems) + 15)
# bars
barplot!(ax_dark, enum_indices, pct_dark, color=:lightblue)
# manual error bars with caps
cap_width = 0.15
for (i, sem) in zip(enum_indices, sems)
    low = pct_dark[i] - sem
    high = pct_dark[i] + sem
    # vertical error line
    lines!(ax_dark, [i, i], [low, high], color=:black)
    # cap on top
    lines!(ax_dark, [i - cap_width, i + cap_width], [high, high], color=:black)
end

# annotate significance for paired comparisons
pairs = [(1,2), (3,4), (5,6), (7,8)]
for (i, j) in pairs
    v1, v2 = rep_pcts_list[i], rep_pcts_list[j]
    if !isempty(v1) && !isempty(v2)
        minSamples = min(length(v1),length(v2))
        # independent two-sample Welch t-test manually
        p = pvalue(OneSampleTTest(v1[1:minSamples], v2[1:minSamples])) 
        # y-position above error bars
        y = maximum([pct_dark[i] + sems[i], pct_dark[j] + sems[j]]) + 5
        # horizontal line
        lines!(ax_dark, [i, j], [y, y], color=:black)
        # vertical ticks
        lines!(ax_dark, [i, i], [y-1, y], color=:black)
        lines!(ax_dark, [j, j], [y-1, y], color=:black)
        # significance label
        # print significance to terminal and draw label
        annot_y = y + 2
        # draw significance label via scatter! with text annotation

        text!(ax_dark,[mean([i, j])],[annot_y],text  = significance_level(p),color = :black,align = (:center, :bottom),fontsize = 28)

    end
end

# display
display(GLMakie.Screen(), fig_dark)