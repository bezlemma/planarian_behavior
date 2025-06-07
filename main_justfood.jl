include("Stab_Interp_Track_Worms.jl")
using .Stab_Interp_Track_Worms: track_worms


filepaths_test = [ "/Users/dl0346/Documents/PlanarianVideos/May28/RedFood/ATMUNC4_1.tif", 
    "/Users/dl0346/Documents/PlanarianVideos/May28/RedFood/ATMUNC4_2.tif", 
    "/Users/dl0346/Documents/PlanarianVideos/May28/RedFood/ATMUNC4_3.tif", 
    "/Users/dl0346/Documents/PlanarianVideos/May28/RedFood/ATMUNC4_4.tif", 
    "/Users/dl0346/Documents/PlanarianVideos/May28/RedFood/ATMUNC4_5.tif", 
    "/Users/dl0346/Documents/PlanarianVideos/May28/RedFood/ATMUNC4_6.tif",
]
thresholds_test = [1.05, 1.05, 1.05, 1.05, 1.05, 1.05] # Represents % of brightest pixels to be selected
filepaths = [
    "/Users/dl0346/Documents/PlanarianVideos/May28/RedFood/ATM1_1.tif", 
    "/Users/dl0346/Documents/PlanarianVideos/May28/RedFood/ATM1_2.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/RedFood/ATM1_3.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/RedFood/ATM1_4.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/RedFood/ATM1_5.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/RedFood/ATM1_6.tif",

    "/Users/dl0346/Documents/PlanarianVideos/May28/RedFood/ATMUNC1_1.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/RedFood/ATMUNC1_2.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/RedFood/ATMUNC1_3.tif", 
    "/Users/dl0346/Documents/PlanarianVideos/May28/RedFood/ATMUNC1_4.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/RedFood/ATMUNC1_5.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/RedFood/ATMUNC1_6.tif", 

    "/Users/dl0346/Documents/PlanarianVideos/May28/RedFood/ATM4_1.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/RedFood/ATM4_2.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/RedFood/ATM4_3.tif", 
    "/Users/dl0346/Documents/PlanarianVideos/May28/RedFood/ATM4_4.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/RedFood/ATM4_5.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/RedFood/ATM4_6.tif", 

    "/Users/dl0346/Documents/PlanarianVideos/May28/RedFood/ATMUNC4_1.tif", 
    "/Users/dl0346/Documents/PlanarianVideos/May28/RedFood/ATMUNC4_2.tif", 
    "/Users/dl0346/Documents/PlanarianVideos/May28/RedFood/ATMUNC4_3.tif", 
    "/Users/dl0346/Documents/PlanarianVideos/May28/RedFood/ATMUNC4_4.tif", 
    "/Users/dl0346/Documents/PlanarianVideos/May28/RedFood/ATMUNC4_5.tif", 
    "/Users/dl0346/Documents/PlanarianVideos/May28/RedFood/ATMUNC4_6.tif", 

    "/Users/dl0346/Documents/PlanarianVideos/May28/RedFood/ATM6_1.tif", 
    "/Users/dl0346/Documents/PlanarianVideos/May28/RedFood/ATM6_2.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/RedFood/ATM6_3.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/RedFood/ATM6_4.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/RedFood/ATM6_5.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/RedFood/ATM6_6.tif", 

    "/Users/dl0346/Documents/PlanarianVideos/May28/RedFood/NDK1_1.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/RedFood/NDK1_2.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/RedFood/NDK1_3.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/RedFood/NDK1_4.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/RedFood/NDK1_5.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/RedFood/NDK1_6.tif",

    "/Users/dl0346/Documents/PlanarianVideos/May28/RedFood/NDKUNC1_1.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/RedFood/NDKUNC1_2.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/RedFood/NDKUNC1_3.tif", 
    "/Users/dl0346/Documents/PlanarianVideos/May28/RedFood/NDKUNC1_4.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/RedFood/NDKUNC1_5.tif", 
    "/Users/dl0346/Documents/PlanarianVideos/May28/RedFood/NDKUNC1_6.tif", 

    "/Users/dl0346/Documents/PlanarianVideos/May28/RedFood/TEC1_1.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/RedFood/TEC1_2.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/RedFood/TEC1_3.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/RedFood/TEC1_4.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/RedFood/TEC1_5.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/RedFood/TEC1_6.tif",

    "/Users/dl0346/Documents/PlanarianVideos/May28/RedFood/TECUNC1_1.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/RedFood/TECUNC1_2.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/RedFood/TECUNC1_3.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/RedFood/TECUNC1_4.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/RedFood/TECUNC1_5.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/RedFood/TECUNC1_6.tif",
]

thresholds = [0.05, 0.05, 0.05, 0.05, 0.05, 0.05,  #ATM1
              0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 
              0.05, 0.05, 0.05, 0.05, 0.05, 0.05,  #ATM4
              1.05, 1.05, 1.05, 1.05, 1.05, 1.05,  #ATMUNC4
            0.05, 0.05, 0.05, 0.05, 0.05, 0.05,    #ATM6
            0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 
            0.05, 0.05, 0.05, 0.05, 0.05, 0.05,    #NDKUNC
            0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 
            0.05, 0.05, 0.05, 0.05, 0.05, 0.05] 

results, stacks = track_worms(filepaths; pixel_percentages=thresholds)


include("Plotting_Aux.jl")
using .Plotting_Aux
#plot_lightdark(results) # Generate light/dark plot
view_stack_and_worm(stacks[end],results[end]) #Use the live-viewer

##

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
conditions = ["ATM 1wpi", "Control", "ATM 4wpi", "Control", "ATM 6wpi","NDK 1wpi", "Control", "TEC 1wpi", "Control"]
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
        ys = [pt[1] for pt in results[idx].positions_cm]
        thr = results[idx].img_rows * 1.0 / 58.6 / 2 #THIS CONVERTS TO PIXEL/CM AND HTEN IS HALF THAT! THIS SOHULDN"T BE HANDLED HERE!!!!
        pct = count(y -> y < thr, ys) / length(ys) * 100
        push!(rep_pcts, pct)
    end
    # store replicate percentages
    push!(rep_pcts_list, rep_pcts)
    # compute and store mean
    push!(pct_dark, mean(rep_pcts))
end

# Plot bar graph with error bars and significance
function significance_level(p)
    if p < 0.001; println("p=$p, ***"); return "***"
    elseif p < 0.01; println("p=$p, **"); return "**"
    elseif p < 0.05; println("p=$p, *"); return "*"
    else println("p=$p, ns"); return "ns"; end
end

# compute standard deviation of the mean for each condition
stds = [ isempty(v) ? NaN : std(v) for v in rep_pcts_list ]

fig_dark = Figure()
ax_dark = Axis(fig_dark[1,1], xticks=(enum_indices, conditions), ylabel="% Time in Non-Food Half", xticklabelrotation=45)

# set y limits to cover error bars and room for annotations
ylims!(ax_dark, 0, maximum(pct_dark .+ stds) + 25) # Increased spacing for annotations

# bar plot
barplot!(ax_dark, enum_indices, pct_dark, 
color=[:dodgerblue, :darkblue, :dodgerblue, :darkblue, :dodgerblue, :orange, :darkorange, :orange, :darkorange])

# error bars
cap_width = 0.15
for (i, sem) in zip(enum_indices, stds)
    low = pct_dark[i] - sem
    high = pct_dark[i] + sem
    lines!(ax_dark, [i, i], [low, high], color=:black)
end


pairs = [(3,4)]
occupied_levels = [] # Stores (x_start, x_end, y_pos) of existing bars
vertical_step = 0.5  # Vertical distance to shift bars to avoid overlap

#We have a lot of logic here to try to place significance values in locations that don't overlap with each other
for (i, j) in pairs
    v1, v2 = rep_pcts_list[i], rep_pcts_list[j]
    if !isempty(v1) && !isempty(v2)
        minSamples = min(length(v1),length(v2))
        # independent two-sample Welch t-test
        p = pvalue(OneSampleTTest(v1[1:minSamples], v2[1:minSamples]))

        y = maximum([pct_dark[i] + stds[i], pct_dark[j] + stds[j]]) + 5
        
        max_iterations = 4 # Safety break to prevent any possibility of an infinite loop
        for iter in 1:max_iterations
            max_conflicting_y = -Inf
            for (x_start, x_end, y_level) in occupied_levels
                x_overlaps = max(i, j) > x_start && min(i, j) < x_end
                y_conflicts = abs(y - y_level) < vertical_step
                
                if x_overlaps && y_conflicts
                    max_conflicting_y = max(max_conflicting_y, y_level)
                end
            end

            if max_conflicting_y > -Inf
                y = max_conflicting_y + vertical_step
            else
                break 
            end
            if iter == max_iterations
                @warn "Could not resolve annotation overlap for pair ($i, $j) after $max_iterations iterations. Placing annotation at last calculated position."
            end
        end
        
        # Add the new bar's final, non-conflicting position to the list
        push!(occupied_levels, (min(i,j), max(i,j), y))
        
        # Draw the significance annotation at the calculated `y`
        lines!(ax_dark, [i, j], [y, y], color=:black)
        lines!(ax_dark, [i, i], [y-1, y], color=:black)
        lines!(ax_dark, [j, j], [y-1, y], color=:black)
        text!(ax_dark,[mean([i, j])],[y], text = significance_level(p), color = :black, align = (:center, :bottom), fontsize = 28)
    end
end

display(GLMakie.Screen(), fig_dark)


