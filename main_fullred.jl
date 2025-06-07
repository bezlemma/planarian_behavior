include("Stab_Interp_Track_Worms.jl")
using .Stab_Interp_Track_Worms: track_worms


filepaths_test = ["/Users/dl0346/Documents/PlanarianVideos/May28/FullRed/ATM1_1.tif"]
thresholds_test = [0.05] # Represents % of brightest pixels to be selected
filepaths = [
    "/Users/dl0346/Documents/PlanarianVideos/May28/FullRed/ATM1_1.tif", 
    "/Users/dl0346/Documents/PlanarianVideos/May28/FullRed/ATM1_2.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/FullRed/ATM1_3.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/FullRed/ATM1_4.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/FullRed/ATM1_5.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/FullRed/ATM1_6.tif",

    "/Users/dl0346/Documents/PlanarianVideos/May28/FullRed/ATMUNC1_1.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/FullRed/ATMUNC1_2.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/FullRed/ATMUNC1_3.tif", 
    "/Users/dl0346/Documents/PlanarianVideos/May28/FullRed/ATMUNC1_4.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/FullRed/ATMUNC1_5.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/FullRed/ATMUNC1_6.tif", 

    "/Users/dl0346/Documents/PlanarianVideos/May28/FullRed/ATM4_1.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/FullRed/ATM4_2.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/FullRed/ATM4_3.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/FullRed/ATM4_4.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/FullRed/ATM4_5.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/FullRed/ATM4_6.tif", 

    "/Users/dl0346/Documents/PlanarianVideos/May28/FullRed/ATMUNC4_1.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/FullRed/ATMUNC4_2.tif", 
    "/Users/dl0346/Documents/PlanarianVideos/May28/FullRed/ATMUNC4_3.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/FullRed/ATMUNC4_4.tif", 
    "/Users/dl0346/Documents/PlanarianVideos/May28/FullRed/ATMUNC4_5.tif", 
    "/Users/dl0346/Documents/PlanarianVideos/May28/FullRed/ATMUNC4_6.tif", 
    

    "/Users/dl0346/Documents/PlanarianVideos/May28/FullRed/ATM6_1.tif", 
    "/Users/dl0346/Documents/PlanarianVideos/May28/FullRed/ATM6_2.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/FullRed/ATM6_3.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/FullRed/ATM6_4.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/FullRed/ATM6_5.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/FullRed/ATM6_6.tif", 

    "/Users/dl0346/Documents/PlanarianVideos/May28/FullRed/NDK1_1.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/FullRed/NDK1_2.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/FullRed/NDK1_3.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/FullRed/NDK1_4.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/FullRed/NDK1_5.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/FullRed/NDK1_6.tif",

    "/Users/dl0346/Documents/PlanarianVideos/May28/FullRed/NDKUNC1_1.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/FullRed/NDKUNC1_2.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/FullRed/NDKUNC1_3.tif", 
    "/Users/dl0346/Documents/PlanarianVideos/May28/FullRed/NDKUNC1_4.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/FullRed/NDKUNC1_5.tif", 
    "/Users/dl0346/Documents/PlanarianVideos/May28/FullRed/NDKUNC1_6.tif",

    "/Users/dl0346/Documents/PlanarianVideos/May28/FullRed/TEC1_1.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/FullRed/TEC1_2.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/FullRed/TEC1_3.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/FullRed/TEC1_4.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/FullRed/TEC1_5.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/FullRed/TEC1_6.tif",

    "/Users/dl0346/Documents/PlanarianVideos/May28/FullRed/TECUNC1_1.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/FullRed/TECUNC1_2.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/FullRed/TECUNC1_3.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/FullRed/TECUNC1_4.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/FullRed/TECUNC1_5.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/FullRed/TECUNC1_6.tif",
]

thresholds = [0.05, 0.05, 0.05, 0.05, 0.05, 0.05,  #ATM1
              0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 
              0.05, 0.05, 0.05, 0.05, 0.05, 0.05,  #ATM4
              0.05, 0.05, 0.05, 0.05, 0.05, 0.05,  #ATMUNC4
            0.05, 0.05, 0.05, 0.05, 0.05, 0.05,   #ATM6
            0.05, 0.05, 0.05, 0.05, 0.05, 0.05,
            0.05, 0.05, 0.05, 0.05, 0.05, 0.05,    #NDKUNC
            0.05, 0.05, 0.05, 0.05, 0.05, 0.05,
            0.05, 0.05, 0.05, 0.05, 0.05, 0.05] 

results, stacks = track_worms(filepaths; pixel_percentages=thresholds)


include("Plotting_Aux.jl")
using .Plotting_Aux
#plot_lightdark(results) # Generate light/dark plot
#view_stack_and_worm(stacks[end],results[end]) #Use the live-viewer

cond_labels = [  "ATM 1wpi" , "control", "ATM 4wip", "control", "ATM 6wpi", "NDK 1wpi", "control", "TEC 1wpi", "control"] #Put whatever conditions you would like to plot here
plot_behavior_distribution(results, cond_labels)

##
using GLMakie, Statistics, HypothesisTests, Colors
#
# Add diffusion behavior plot by averaging MSDs per condition
using .Plotting_Aux: calculate_msd
using GeometryBasics: Point2f

fig = Figure(size=(1000, 750))
ax = Axis(fig[1, 1],
    xlabel="Time Lag τ (s)", ylabel="MSD (cm²)",
    xscale=log10, yscale=log10,
    title="Average MSD per Condition")
colors = distinguishable_colors(length(cond_labels), [RGB(1,1,1), RGB(0,0,0)], lchoices=range(20, stop=70, length=15))
global_lags = Float32[]  # track all lags across conditions
first_vals = Float32[]  # track first MSD values for anchor
for (i, cond) in enumerate(cond_labels)
    cond_data = filter(d -> occursin(cond, d.filepath), results)
    if isempty(cond_data)
        continue
    end
    # compute individual MSDs
    msd_list = [calculate_msd([Point2f(p[1], p[2]) for p in d.positions_cm], d.times_s) for d in cond_data]
    # aggregate lags
    all_lags = sort(unique(vcat([collect(keys(m)) for m in msd_list]...)))
    mean_msd = [mean([m[lag] for m in msd_list if haskey(m, lag)]) for lag in all_lags]
    lines!(ax, all_lags, mean_msd, color=colors[i], linewidth=2, label=cond)
    # accumulate for reference lines
    append!(global_lags, all_lags)
    if !isempty(mean_msd)
        push!(first_vals, mean_msd[1])
    end
end
# add reference slopes if we have data
if !isempty(global_lags) && !isempty(first_vals)
    min_lag, max_lag = minimum(global_lags), maximum(global_lags)
    ref_lags = 10 .^ LinRange(log10(min_lag), log10(max_lag), 50)
    anchor = minimum(first_vals)
    lines!(ax, ref_lags, anchor*(ref_lags ./ ref_lags[1]).^1, color=RGB(0.5,0.5,0.5), linestyle=:dash, label="α=1")
    lines!(ax, ref_lags, anchor*(ref_lags ./ ref_lags[1]).^0.5, color=RGB(0.3,0.3,0.3), linestyle=:dashdot, label="α=0.5")
    lines!(ax, ref_lags, anchor*(ref_lags ./ ref_lags[1]).^2, color=RGB(0.7,0.7,0.7), linestyle=:dot, label="α=2")
    axislegend(ax, position=:rb)
end

# Bar plot of mean speed per condition
ax2 = Axis(fig[2, 1],
    xlabel="Condition", ylabel="Mean Speed (cm/s)",
    title="Mean Speed per Condition")
# compute per-track mean speeds per condition
speeds_per_cond = [begin
    cond_data = filter(d -> occursin(cond, d.filepath), results)
    [mean([sqrt((d.positions_cm[j+1][1] - d.positions_cm[j][1])^2 +
                (d.positions_cm[j+1][2] - d.positions_cm[j][2])^2) /
           (d.times_s[j+1] - d.times_s[j])
       for j in 1:(length(d.positions_cm)-1)]) for d in cond_data]
end for cond in cond_labels]
mean_speeds = [mean(s) for s in speeds_per_cond]
# draw bars
barplot!(ax2, 1:length(cond_labels), mean_speeds, color=[:dodgerblue, :darkblue, :dodgerblue, :darkblue, :dodgerblue, :orange, :darkorange, :orange, :darkorange])
# overlay individual points
for (i, speeds) in enumerate(speeds_per_cond)
    scatter!(ax2, fill(i, length(speeds)), speeds, color=:black)
end
ax2.xticks = (1:length(cond_labels), cond_labels)
# render updated figure
display(GLMakie.Screen(), fig)



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
        ys = [pt[2] for pt in results[idx].positions_cm]
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
    if p < 0.001; println("p=$p, ★★★"); return "★★"
    elseif p < 0.01; println("p=$p, ★★"); return "★★"
    elseif p < 0.05; println("p=$p, ★"); return "★"
    else println("p=$p, ns"); return "ns"; end
end

# compute standard deviation of the mean for each condition
stds = [ isempty(v) ? NaN : std(v) for v in rep_pcts_list ]

fig_dark = Figure()
ax_dark = Axis(fig_dark[1,1], xticks=(enum_indices, conditions), ylabel="% Time in Bottom Half", xticklabelrotation=45)

# set y limits to cover error bars and room for annotations
ylims!(ax_dark, 0, maximum(pct_dark .+ stds) + 25) # Increased spacing for annotations

# bar plot
barplot!(ax_dark, enum_indices, pct_dark, 
color=[:dodgerblue, :darkblue, :dodgerblue, :darkblue, :dodgerblue, :orange, :darkorange, :orange, :darkorange])

# error bars
for (i, sem) in zip(enum_indices, stds)
    low = pct_dark[i] - sem
    high = pct_dark[i] + sem
    lines!(ax_dark, [i, i], [low, high], color=:black)
end


pairs = [(1,2), (3,4), (6,7), (8,9), 
         (1,3), (1,5), 
         (6,8)]
occupied_levels = [] # Stores (x_start, x_end, y_pos) of existing bars
vertical_step = 8.5  # Vertical distance to shift bars to avoid overlap

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
        text!(ax_dark,[mean([i, j])],[y], text = significance_level(p), color = :black, align = (:center, :bottom), fontsize = 18)
    end
end

display(GLMakie.Screen(), fig_dark)