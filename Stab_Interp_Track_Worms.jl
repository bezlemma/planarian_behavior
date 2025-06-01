#TODO: Smooth out track using a morphosnakes-like approach
#TODO: Start putting together beahviors such as turns
#TODO: Find "jumping" bug and fix

using TiffImages
using GLMakie
using LinearAlgebra # For norm (distance calculation)
using Statistics    # For mean
using Colors        # For distinguishable_colors, RGB, Gray
using GeometryBasics # For Polygon, Rect, Point2f, Point3f, Circle
using ProgressMeter # For progress bar
using Base.Threads  # For multithreading
using Observables   # For interactive Observables with GLMakie

# Import image stabilization utilities
include("FUNC_ImageProcessing.jl")
using .FUNC_ImageProcessing

include("FUNC_ObjectFind.jl")
using .FUNC_ObjectFind

include("FUNC_WormFinderCIRCLE.jl")
using .FUNC_WormFinderCIRCLE

# --- Constants ---, more constants are hiding in FUNC_WormFinderCIRCLE right now
const THRESHOLD = 1.5f0
const CM_PER_PIXEL = 1.0f0 / 58.6f0
const MSEC_PER_FRAME = 50.0f0 * 2.0f0


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

    println("Processing: $filepath")
    raw_data = TiffImages.load(filepath)

    img_rows, img_cols, num_frames = size(raw_data, 1), size(raw_data, 2), size(raw_data, 3)

    data = 1.0f0 .- Float32.(raw_data)
    reference = median(data, dims=3)

    #Comment or uncomment this to stabilize the data.
    #data,reference = FUNC_ImageProcessing.stabilize(data,reference)
    
    divided_stack = data ./ (reference .+ eps(Float32(1.0)))
    global binary_stack = divided_stack  .> THRESHOLD

    positions_, areas_2, major_axes_, minor_axes_, track_img_rows, track_img_cols = create_track(binary_stack)

    times_s_vec = [p[3] * MSEC_PER_FRAME / 1000.0f0 for p in positions_]
    positions_cm_s_vec = [Point3f(p[1]*CM_PER_PIXEL, (track_img_rows - p[2])*CM_PER_PIXEL, times_s_vec[i]) for (i,p) in enumerate(positions_)]

    #Save the WormData struct
    push!(all_worms, WormData(filepath, positions_, areas_2, major_axes_, minor_axes_, img_rows, img_cols, times_s_vec, positions_cm_s_vec))
end


## --------- Plotting below 

## Fun 3D Plot
fig_3d_traj = Figure(size=(1200, 900)) 
max_time_s_overall = 0.0f0; img_width_cm_overall = 7.5f0; img_height_cm_overall = 7.5f0
valid_times_all_worms = [data.times_s for data in all_worms if !isempty(data.times_s)]

max_times_per_dataset = [maximum(ts) for ts in valid_times_all_worms if !isempty(ts)]
first_data_for_dims = all_worms[1]
img_width_cm_overall = first_data_for_dims.img_cols * CM_PER_PIXEL
img_height_cm_overall = first_data_for_dims.img_rows * CM_PER_PIXEL

ax_3d = Axis3(fig_3d_traj[1,1], xlabel="X (cm)", ylabel="Y (cm)", zlabel="Time (s)", 
limits = ( (0, img_width_cm_overall),(0, img_height_cm_overall),(0, maximum(max_times_per_dataset))  ))

colors_traj = distinguishable_colors(max(1, length(all_worms)), [RGB(0.0, 1.0, 1.0), RGB(0.125, 0.875, 0.875), RGB(0.25, 0.75, 0.75), RGB(0.375, 0.625, 0.625), RGB(0.5, 0.5, 0.5)], lchoices=range(20, stop=70, length=15))
for (i, worm) in enumerate(all_worms)
    lines!(ax_3d, worm.positions_cm, color=colors_traj[i], linewidth=4)
end
display(GLMakie.Screen(), fig_3d_traj) 

#interactive_fig = view_stack_and_worm(binary_stack, all_worms[end] ) #if you want to use this, you'll need to load the plot function from Plotting_Aux.jl

