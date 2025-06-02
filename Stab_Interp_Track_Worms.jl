#TODO: Smooth out track using a morphosnakes-like approach
#TODO: Start putting together behaviors such as turns

using TiffImages #Loading
using GLMakie, Observables #Plotting
using LinearAlgebra, Statistics, Colors, GeometryBasics #Basics
using ProgressMeter, Base.Threads  # Interface
using ImageMorphology

# Import image stabilization utilities
include("FUNC_ImageProcessing.jl"); using .FUNC_ImageProcessing

#Import the worm finder script
include("FUNC_WormFinder.jl"); using .FUNC_WormFinder
include("FUNC_ObjectFind.jl"); using .FUNC_ObjectFind
include("FUNC_WormBehaviors.jl"); using .FUNC_WormBehaviors

# --- Constants ---, more constants are hiding in FUNC_WormFinderCIRCLE right now
THRESHOLD = 1.5
CM_PER_PIXEL = 1.0 / 58.6
MSEC_PER_FRAME = 50.0 * 2.0

MIN_AREA = 15 #Min # pixels that are still considered a worm
MAX_AREA = 500 #Max # pixels that are still considered a worm
MAJOR_MINOR_MAX_RATIO = 10.0 #If the worm is 30x longer than it is wide, it isn't a worm.
WORM_DISTANCE_SEARCHADD = 20.0

filepaths = [
    "/Users/dl0346/Documents/PlanarianVideos/May28/TrailFollow/ATM1_1.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/TrailFollow/ATM1_2.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/TrailFollow/ATM1_3.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/TrailFollow/ATM1_4.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/TrailFollow/ATM1_5.tif",
]
global binary_stack = BitArray{3}(undef, (0,0,0))
all_worms = WormData[]

for filepath in filepaths

    #Load the data
    println("Processing: $filepath")
    raw_data = TiffImages.load(filepath)
    img_rows, img_cols, num_frames = size(raw_data, 1), size(raw_data, 2), size(raw_data, 3)

    #Invert the data, and take the average to look at the background
    data = 1.0 .- Float32.(raw_data)
    reference = median(data, dims=3)

    #Comment or uncomment this to stabilize the data.
    #data,reference = FUNC_ImageProcessing.stabilize(data,reference)
    
    #Background divide
    divided_stack = data ./ (reference .+ eps(Float32(1.0)))

    #Segment the data based on a threshold
    binary_mask = divided_stack  .> THRESHOLD

    binary_mask = dilate(erode(binary_mask)) #Remove 1 pixel lines

    #This is enterily redundent, but I find it helpful for debugging to remove too small or too large objects before the full pipeline.
    #This is also cast as a global for easier debugging
    global binary_stack = remove_outlier_objects(binary_mask, MIN_AREA, MAX_AREA, MAJOR_MINOR_MAX_RATIO)

    #Create the track
    circle_status = true
    positions, areas_2, major_axes, minor_axes = create_track(binary_stack, circle_status, WORM_DISTANCE_SEARCHADD)

    #Track in some behaviors, specifically
    behaviors = compute_behaviors(positions, 
                                 [p[3] * MSEC_PER_FRAME / 1000.0 for p in positions],
                                 major_axes, minor_axes,
                                 Point2f(img_cols/2, img_rows/2),
                                 min(img_rows, img_cols)/2)

    times_s = [p[3] * MSEC_PER_FRAME / 1000.0 for p in positions]
    pos_cm = [Point3f(p[1]*CM_PER_PIXEL, (img_rows - p[2])*CM_PER_PIXEL, times_s[i]) for (i,p) in enumerate(positions)]

    #Save the WormData struct (now with behaviors)
    push!(all_worms, WormData(filepath, positions, areas_2, major_axes, minor_axes, img_rows, img_cols, times_s, pos_cm, behaviors))
    
    # Debug: Print behavior distribution
    behavior_counts = Dict{Symbol, Int}()
    for b in behaviors
        behavior_counts[b] = get(behavior_counts, b, 0) + 1
    end
    println("Behavior distribution for $filepath:")
    for (behavior, count) in sort(collect(behavior_counts), by=x->x[2], rev=true)
        percentage = round(count / length(behaviors) * 100, digits=1)
        println("  $behavior: $count ($percentage%)")
    end
end


## --------- Beautiful 3D Trajectory Plotting

# Create a beautiful figure with better styling
fig_3d_traj = Figure(size=(1400, 1000), backgroundcolor=:white) 

valid_times_all_worms = [data.times_s for data in all_worms if !isempty(data.times_s)]
max_time = maximum([maximum(ts) for ts in valid_times_all_worms])
img_width_cm = all_worms[1].img_cols * CM_PER_PIXEL
img_height_cm= all_worms[1].img_rows * CM_PER_PIXEL

# Create 3D axis with better styling
ax_3d = Axis3(fig_3d_traj[1,1], 
    xlabel="X Position (cm)", ylabel="Y Position (cm)", zlabel="Time (s)",
    xlabelsize=14, ylabelsize=14, zlabelsize=14,
    title="Planarian Trajectory with Behavior Classification", titlesize=16,
    limits = ((0, img_width_cm), (0, img_height_cm), (0, max_time)),
    aspect = (1, 1, 0.8)  # Make time axis slightly compressed
)

# Define beautiful colors for each behavior
behavior_colors = Dict(
    :away    => RGBA(1.0, 0.2, 0.2, 0.7),  # Bright Red
    :toward  => RGBA(0.2, 0.5, 1.0, 0.7),  # Bright Blue  
    :along   => RGBA(0.1, 0.1, 0.1, 0.7),  # Black
    :turning => RGBA(0.0, 0.8, 0.2, 0.7),  # Bright Green
    :pausing => RGBA(0.9, 0.1, 0.9, 0.7)   # Bright Magenta
)

# Plot each worm's trajectory as a single smooth continuous line
for (i, worm) in enumerate(all_worms)
    points = worm.positions_cm
    if length(points) < 2
        continue
    end

    # Create smooth color array for the trajectory
    colors = Vector{RGBA{Float32}}(undef, length(points))
    for j in 1:length(points)
        beh = worm.behaviors[j]
        colors[j] = get(behavior_colors, beh, RGBA(0.5,0.5,0.5,0.7)) # Gray default
    end

    # Plot as single thick continuous line with smooth color transitions
    lines!(ax_3d, points, color=colors, linewidth=8.0, linestyle=:solid)
end

# Create a clean legend with sample lines
legend_elements = []
for (beh, color) in behavior_colors
    # Create visible legend lines (not off-screen)
    sample_points = [Point3f(0, 0, max_time * 1.1), Point3f(img_width_cm * 0.1, 0, max_time * 1.1)]
    lines!(ax_3d, sample_points, color=color, linewidth=6.0, visible=false) # Hidden but in legend
    push!(legend_elements, LineElement(color=color, linewidth=4, linestyle=:solid))
end

# Add legend
Legend(fig_3d_traj[1,2], legend_elements, [string(beh) for beh in keys(behavior_colors)], 
       "Behaviors", framevisible=true, backgroundcolor=:white)

display(GLMakie.Screen(), fig_3d_traj) 

#interactive_fig = view_stack_and_worm(binary_stack, all_worms[end] ) #if you want to use this, you'll need to load the plot function from Plotting_Aux.jl