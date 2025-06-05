include("Stab_Interp_Track_Worms.jl")
using .Stab_Interp_Track_Worms: track_worms


filepaths_test = ["/Users/dl0346/Documents/PlanarianVideos/May28/HalfRedFood/NDKUNC1_3.tif"]
thresholds_test = [0.25] # Represents 0.05% of brightest pixels to be selected
filepaths = [
    "/Users/dl0346/Documents/PlanarianVideos/May28/HalfRedFood/ATM1_1.tif", 
    "/Users/dl0346/Documents/PlanarianVideos/May28/HalfRedFood/ATM1_2.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/HalfRedFood/ATM1_3.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/HalfRedFood/ATM1_4.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/HalfRedFood/ATM1_5.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/HalfRedFood/ATM1_6.tif", #Good at 0.

    "/Users/dl0346/Documents/PlanarianVideos/May28/HalfRedFood/ATMUNC1_1.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/HalfRedFood/ATMUNC1_2.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/HalfRedFood/ATMUNC1_3.tif", 
    "/Users/dl0346/Documents/PlanarianVideos/May28/HalfRedFood/ATMUNC1_4.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/HalfRedFood/ATMUNC1_5.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/HalfRedFood/ATMUNC1_6.tif",

    "/Users/dl0346/Documents/PlanarianVideos/May28/HalfRedFood/ATM4_1.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/HalfRedFood/ATM4_2.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/HalfRedFood/ATM4_3.tif", #Good, tested
    "/Users/dl0346/Documents/PlanarianVideos/May28/HalfRedFood/ATM4_4.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/HalfRedFood/ATM4_5.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/HalfRedFood/ATM4_6.tif",

    "/Users/dl0346/Documents/PlanarianVideos/May28/HalfRedFood/ATMUNC4_1.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/HalfRedFood/ATMUNC4_2.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/HalfRedFood/ATMUNC4_3.tif", 
    "/Users/dl0346/Documents/PlanarianVideos/May28/HalfRedFood/ATMUNC4_4.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/HalfRedFood/ATMUNC4_5.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/HalfRedFood/ATMUNC4_6.tif",

    "/Users/dl0346/Documents/PlanarianVideos/May28/HalfRedFood/ATM6_1.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/HalfRedFood/ATM6_2.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/HalfRedFood/ATM6_3.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/HalfRedFood/ATM6_4.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/HalfRedFood/ATM6_5.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/HalfRedFood/ATM6_6.tif",

    "/Users/dl0346/Documents/PlanarianVideos/May28/HalfRedFood/NDK1_1.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/HalfRedFood/NDK1_2.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/HalfRedFood/NDK1_3.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/HalfRedFood/NDK1_4.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/HalfRedFood/NDK1_5.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/HalfRedFood/NDK1_6.tif",

    "/Users/dl0346/Documents/PlanarianVideos/May28/HalfRedFood/NDKUNC1_1.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/HalfRedFood/NDKUNC1_2.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/HalfRedFood/NDKUNC1_3.tif", #Good at 0.25
    "/Users/dl0346/Documents/PlanarianVideos/May28/HalfRedFood/NDKUNC1_4.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/HalfRedFood/NDKUNC1_5.tif", #Good at 0.05
    "/Users/dl0346/Documents/PlanarianVideos/May28/HalfRedFood/NDKUNC1_6.tif", #Good at 0.25

    "/Users/dl0346/Documents/PlanarianVideos/May28/HalfRedFood/TEC1_1.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/HalfRedFood/TEC1_2.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/HalfRedFood/TEC1_3.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/HalfRedFood/TEC1_4.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/HalfRedFood/TEC1_5.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/HalfRedFood/TEC1_6.tif",

    "/Users/dl0346/Documents/PlanarianVideos/May28/HalfRedFood/TECUNC1_1.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/HalfRedFood/TECUNC1_2.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/HalfRedFood/TECUNC1_3.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/HalfRedFood/TECUNC1_4.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/HalfRedFood/TECUNC1_5.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/HalfRedFood/TECUNC1_6.tif",
]

thresholds = [0.05, 0.05, 0.05, 0.05, 0.25, 0.05,  
              0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 
              0.05, 0.05, 0.05, 0.25, 0.05, 0.05, 
              0.25, 0.05, 0.05, 0.25, 0.05, 0.05,  
 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 
            0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 
            0.05, 0.05, 0.25, 0.05, 0.05, 0.25, #NDKUNC
            0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 
            0.25, 0.05, 0.05, 0.05, 0.25, 0.05] 

results, stacks = track_worms(filepaths; pixel_percentages=thresholds)



include("Plotting_Aux.jl")
using .Plotting_Aux
#plot_lightdark(results) # Generate light/dark plot
view_stack_and_worm(stacks[end],results[end]) #Use the live-viewer

using GLMakie, Statistics, HypothesisTests, Colors

## The following is helpful for debugging what is x and what is y 
#fig_dark = Figure()
#ax = Axis(fig_dark[1,1], ylabel="y pos",xlabel="time")
#y_coord = [pt[2] for pt in results[end].positions_cm]
#x_coord = [pt[1] for pt in results[end].positions_cm]
#lines!(ax, 1:length(x_coord),x_coord)
#fig_dark
##


# Compute and plot % time spent in dark (lower half) per condition
# Define conditions and grouping
conditions = ["ATM1", "ATMUNC1", "ATM4", "ATMUNC4", "ATM6","NDK1", "NDKUNC1", "TEC1", "TECUNC1"]
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
ylims!(ax_dark, 0, maximum(pct_dark .+ sems) + 25) # Increased spacing for annotations

# bar plot
barplot!(ax_dark, enum_indices, pct_dark, 
color=[:dodgerblue, :darkblue, :dodgerblue, :darkblue, :dodgerblue, :orange, :darkorange, :orange, :darkorange])

# error bars
cap_width = 0.15
for (i, sem) in zip(enum_indices, sems)
    low = pct_dark[i] - sem
    high = pct_dark[i] + sem
    lines!(ax_dark, [i, i], [low, high], color=:black)
end

pairs = [(1,2), (3,4), (6,7), (8,9), (1,3), (1,5), (6,8)]
occupied_levels = [] # Stores (x_start, x_end, y_pos) of existing bars
vertical_step = 11  # Vertical distance to shift bars to avoid overlap

for (i, j) in pairs
    v1, v2 = rep_pcts_list[i], rep_pcts_list[j]
    if !isempty(v1) && !isempty(v2)
        minSamples = min(length(v1),length(v2))
        # independent two-sample Welch t-test
        p = pvalue(OneSampleTTest(v1[1:minSamples], v2[1:minSamples]))

        # Initial y-position above the highest error bar in the pair
        y = maximum([pct_dark[i] + sems[i], pct_dark[j] + sems[j]]) + 5
        
        # Check against existing significance bars and adjust `y` if needed
        while true
            overlapping = false
            for (x_start, x_end, y_level) in occupied_levels
                # Check if the x-range of the new bar overlaps with an existing one
                if max(i, j) > x_start && min(i, j) < x_end
                    # Check if y-positions are too close
                    if abs(y - y_level) < vertical_step
                        y = y_level + vertical_step # Shift new bar up
                        overlapping = true
                        break # Restart check with the new `y`
                    end
                end
            end
            if !overlapping
                break # Found a clear spot
            end
        end
        
        # Add the new bar's position to the list of occupied levels
        push!(occupied_levels, (min(i,j), max(i,j), y))
        
        # Draw the significance annotation at the calculated `y`
        lines!(ax_dark, [i, j], [y, y], color=:black)
        lines!(ax_dark, [i, i], [y-1, y], color=:black)
        lines!(ax_dark, [j, j], [y-1, y], color=:black)
        text!(ax_dark,[mean([i, j])],[y], text = significance_level(p), color = :black, align = (:center, :bottom), fontsize = 28)
    end
end

display(GLMakie.Screen(), fig_dark)


