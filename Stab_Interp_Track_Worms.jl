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
  #  "/Users/dl0346/Documents/PlanarianVideos/May28/TrailFollow/ATM1_2.tif",
  #  "/Users/dl0346/Documents/PlanarianVideos/May28/TrailFollow/ATM1_3.tif",
  #  "/Users/dl0346/Documents/PlanarianVideos/May28/TrailFollow/ATM1_4.tif",
  #  "/Users/dl0346/Documents/PlanarianVideos/May28/TrailFollow/ATM1_5.tif",
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
end


## --------- Plotting below 

## Fun 3D Plot
fig_3d_traj = Figure(size=(1200, 900)) 

valid_times_all_worms = [data.times_s for data in all_worms if !isempty(data.times_s)]
max_time = maximum([maximum(ts) for ts in valid_times_all_worms])
img_width_cm = all_worms[1].img_cols * CM_PER_PIXEL
img_height_cm= all_worms[1].img_rows * CM_PER_PIXEL

ax_3d = Axis3(fig_3d_traj[1,1], xlabel="X (cm)", ylabel="Y (cm)", zlabel="Time (s)", 
limits = ( (0, img_width_cm),(0, img_height_cm),(0, max_time)  ))

# Define attributes for each behavior (color, linewidth)
behavior_attributes = Dict(
    :away    => (color=RGBA(1.0, 0.4, 0.4, 0.8), linewidth=3.0),  # Lighter Red
    :toward  => (color=RGBA(0.4, 0.4, 1.0, 0.8), linewidth=3.0),  # Lighter Blue
    :along   => (color=RGBA(0.2, 0.2, 0.2, 1.0), linewidth=5.0),  # Dark Gray
    :turning => (color=RGBA(0.0, 0.8, 0.0, 1.0), linewidth=7.0),  # Bright Green
    :pausing => (color=RGBA(0.8, 0.0, 0.8, 1.0), linewidth=7.0)   # Bright Magenta
)

# Plot each worm's trajectory as a single continuous line with varying attributes
for (i, worm) in enumerate(all_worms)
    points = worm.positions_cm
    if length(points) < 2
        continue # Not enough points to draw a line
    end

    point_colors = Vector{RGBA{Float32}}(undef, length(points))
    point_linewidths = Vector{Float32}(undef, length(points))

    for j in 1:length(points)
        beh = worm.behaviors[j]
        attrs = get(behavior_attributes, beh, (color=RGBA(0.5,0.5,0.5,0.5), linewidth=2.0)) # Default for unknown
        point_colors[j] = attrs.color
        point_linewidths[j] = attrs.linewidth
    end

    lines!(ax_3d, points, color=point_colors, linewidth=point_linewidths, linestyle=:solid)
end

# Legend: dummy lines for each behavior style
for (beh, attrs) in behavior_attributes
    dummy = [Point3f(-1,-1,-1), Point3f(-2,-2,-2)] # Off-screen points
    lines!(ax_3d, dummy, color=attrs.color, linewidth=attrs.linewidth, linestyle=:solid, label=string(beh))
end
axislegend(ax_3d, position=:rb)

display(GLMakie.Screen(), fig_3d_traj) 

#interactive_fig = view_stack_and_worm(binary_stack, all_worms[end] ) #if you want to use this, you'll need to load the plot function from Plotting_Aux.jl

##
# Ensure all_worms, CM_PER_PIXEL, Point3f, RGBA, Figure, Axis3,
# lines!, Legend, LineElement, display, GLMakie.Screen,
# and variables like img_width_cm, img_height_cm, max_time are defined in your scope.

fig_3d_traj = Figure(size=(1200, 900))

# Assuming img_width_cm, img_height_cm, max_time are available from your setup
ax_3d = Axis3(fig_3d_traj[1,1], xlabel="X (cm)", ylabel="Y (cm)", zlabel="Time (s)",
              limits = ( (0, img_width_cm),(0, img_height_cm),(0, max_time) ))

# --- Behavior Attributes and Colormap Setup ---
behavior_attributes = Dict(
    :away    => (color=RGBA(1.0f0, 0.4f0, 0.4f0, 0.8f0), linewidth=7.0f0),
    :toward  => (color=RGBA(0.4f0, 0.4f0, 1.0f0, 0.8f0), linewidth=7.0f0),
    :along   => (color=RGBA(0.2f0, 0.2f0, 0.2f0, 1.0f0), linewidth=7.0f0),
    :turning => (color=RGBA(0.0f0, 0.8f0, 0.0f0, 1.0f0), linewidth=7.0f0),
    :pausing => (color=RGBA(0.8f0, 0.0f0, 0.8f0, 1.0f0), linewidth=7.0f0)
)
default_attrs = (color=RGBA(0.2f0,0.2f0,0.2f0,0.2f0), linewidth=7.0f0)
# Symbol for behaviors not in behavior_attributes or for mismatched data
default_category_symbol = :default_or_unknown_behavior 

# Define an explicit order for behaviors to ensure consistent ID mapping for the colormap
ordered_behavior_keys = [:away, :toward, :along, :turning, :pausing]

# Build the map from behavior symbol to a 1-based integer ID
# Also, build the colormap array (vector of colors) based on this order
behavior_to_id = Dict{Symbol, Int}()
colormap_for_lines = RGBA{Float32}[] # This will be our custom colormap

current_id = 1
for beh_key in ordered_behavior_keys
    if haskey(behavior_attributes, beh_key)
        behavior_to_id[beh_key] = current_id
        push!(colormap_for_lines, behavior_attributes[beh_key].color)
        current_id += 1
    # else: if a key in ordered_behavior_keys is not in behavior_attributes, it's skipped
    end
end

# Add a category and color for any default/unknown behaviors
behavior_to_id[default_category_symbol] = current_id
push!(colormap_for_lines, default_attrs.color)
num_total_categories = current_id # Total number of distinct categories for colorrange

# --- Plotting Trajectories ---
for worm in all_worms # Assuming all_worms is an iterable of your worm data structures
    points = worm.positions_cm
    behaviors = worm.behaviors # Assuming each worm object has these fields

    if length(points) < 2
        continue
    end

    point_behavior_ids = Vector{Int}(undef, length(points))
    point_linewidths = Vector{Float32}(undef, length(points))

    for j in 1:length(points)
        # Determine the behavior symbol for the current point
        # If behaviors array is shorter than points, or behavior is not recognized, use default
        current_behavior_symbol = (j <= length(behaviors)) ? behaviors[j] : default_category_symbol
        
        # Get the integer ID for this behavior, defaulting to the ID for unknown/default
        point_behavior_ids[j] = get(behavior_to_id, current_behavior_symbol, behavior_to_id[default_category_symbol])
        
        # Determine linewidth based on the original behavior attributes
        attrs_for_linewidth = get(behavior_attributes, current_behavior_symbol, default_attrs)
        point_linewidths[j] = attrs_for_linewidth.linewidth
    end

    # Plot using the integer IDs for color, and provide the custom colormap and colorrange
    lines!(ax_3d, points,
           color = point_behavior_ids,          # Array of integer IDs
           colormap = colormap_for_lines,       # Custom array of RGBA colors
           colorrange = (1, num_total_categories), # Maps IDs 1..N to the colormap
           linewidth = point_linewidths,
           linestyle = :solid,
           ssao = true) # Optional: Screen Space Ambient Occlusion for better depth perception
end

# --- Legend ---
legend_elements = Makie.LegendElement[] # Specify type for clarity
legend_labels = String[]

# Create legend entries for the defined behaviors
for beh_key in ordered_behavior_keys
    if haskey(behavior_to_id, beh_key) && haskey(behavior_attributes, beh_key)
        id = behavior_to_id[beh_key]
        # Ensure the ID is within the bounds of the colormap
        if 1 <= id <= length(colormap_for_lines)
            actual_color_for_legend = colormap_for_lines[id]
            attrs_for_legend = behavior_attributes[beh_key] # For consistent linewidth in legend
            
            push!(legend_elements, LineElement(color = actual_color_for_legend, linestyle = :solid, linewidth = attrs_for_legend.linewidth))
            push!(legend_labels, string(beh_key))
        end
    end
end

# Optionally, add a legend entry for the default/unknown category if desired
# id_default = behavior_to_id[default_category_symbol]
# if 1 <= id_default <= length(colormap_for_lines)
#     push!(legend_elements, LineElement(color = colormap_for_lines[id_default], linestyle = :solid, linewidth = default_attrs.linewidth))
#     push!(legend_labels, string(default_category_symbol))
# end

# Add legend to the axis, similar to your original approach
if !isempty(legend_elements)
    axislegend(ax_3d, legend_elements, legend_labels, "Behavior", position = :rb)
end

display(GLMakie.Screen(), fig_3d_traj)