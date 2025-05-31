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

# --- Constants ---
const THRESHOLD = 1.5f0
const MIN_AREA = 15
const WORM_DISTANCE_LINK = 40.0f0       
const WORM_DISTANCE_SEARCHADD = 0.1f0  
const MIN_ACCEPTABLE_CIRCULARITY = 0.10f0
const MAX_AREA = 500 \

const CM_PER_PIXEL = 1.0f0 / 58.6f0
const MSEC_PER_FRAME = 50.0f0 * 2.0f0

# --- Helper Functions (Original) ---
# tracking functions and helpers now imported from FUNC_WormFinderCIRCLE module
function view_binary_stack_with_tracks_glmakie(
    binary_stack_to_show::BitArray{3},
    tracked_worm_data::WormData;
    marker_color=RGBAf(1.0, 0.0, 0.0, 0.7) 
)
    rows, cols, num_frames = size(binary_stack_to_show)
    worm_info_for_frame = Dict{Int, Tuple{Point2f, Int}}() 
    for i in 1:length(tracked_worm_data.positions_)
        pt3 = tracked_worm_data.positions_[i]
        frame_idx = Int(round(pt3[3]))
        centroid = Point2f(pt3[2], pt3[1]) # (col_, row_) for GLMakie point
        area = tracked_worm_data.areas_2[i]
        worm_info_for_frame[frame_idx] = (centroid, area)
    end
    fig = Figure(size = (cols > rows ? (800, 800 * rows/cols + 100) : (800 * cols/rows + 100, 800)))
    ax_img = Axis(fig[1, 1])
    slider = Slider(fig[2, 1], range=1:num_frames, startvalue=1)
    frame_idx_obs = slider.value
    img_slice_obs = lift(frame_idx_obs) do f_idx
        return Gray.(view(binary_stack_to_show, :, :, f_idx)) 
    end
    image!(ax_img, img_slice_obs, interpolate=false, colormap=:grays, colorrange=(0,1))
    worm_marker_obs = Observable([Circle(Point2f(NaN, NaN), 0f0)])
    on(frame_idx_obs) do f_idx
        if haskey(worm_info_for_frame, f_idx)
            centroid, wormarea = worm_info_for_frame[f_idx]
            radius = sqrt(max(0.0, Float64(wormarea)) / π)
            display_radius = clamp(Float32(radius), 2f0, 25f0)
            worm_marker_obs[] = [Circle(centroid, display_radius)]
        else
            worm_marker_obs[] = [Circle(Point2f(NaN, NaN), 0f0)]
        end
    end
    poly!(ax_img, worm_marker_obs, color=marker_color, strokecolor=:transparent)
    display(GLMakie.Screen(), fig)
    return fig
end


filepaths = [
  #   "/Users/dl0346/Documents/PlanarianVideos/May28/TrailFollow/ATM1_1.tif",
    "/Users/dl0346/Documents/PlanarianVideos/May28/TrailFollow/ATM1_2.tif",
  # "/Users/dl0346/Documents/PlanarianVideos/May28/TrailFollow/ATM1_3.tif",
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

    positions_, areas_2, perimeters_, circularities, track_img_rows, track_img_cols = create_track(binary_stack)

    times_s_vec = [p[3] * MSEC_PER_FRAME / 1000.0f0 for p in positions_]
    positions_cm_s_vec = [Point3f(p[1]*CM_PER_PIXEL, (track_img_rows - p[2])*CM_PER_PIXEL, times_s_vec[i]) for (i,p) in enumerate(positions_)]

    #Save the WormData struct
    push!(all_worms, WormData(filepath, positions_, areas_2, perimeters_, circularities, img_rows, img_cols, times_s_vec, positions_cm_s_vec))
end

## 3D Plot
fig_3d_traj = Figure(size=(1200, 900)) 
max_time_s_overall = 0.0f0; img_width_cm_overall = 7.5f0; img_height_cm_overall = 7.5f0
valid_times_all_worms = [data.times_s for data in all_worms if !isempty(data.times_s)]

max_times_per_dataset = [maximum(ts) for ts in valid_times_all_worms if !isempty(ts)]
first_data_for_dims = all_worms[1]
img_width_cm_overall = first_data_for_dims.img_cols * CM_PER_PIXEL
img_height_cm_overall = first_data_for_dims.img_rows * CM_PER_PIXEL

ax_3d = Axis3(fig_3d_traj[1,1], xlabel="X (cm)", ylabel="Y (cm)", zlabel="Time (s)", 
limits = ( (0, img_width_cm_overall),(0, img_height_cm_overall),(0, maximum(max_times_per_dataset))  ))

colors_traj = distinguishable_colors(max(1,length(all_worms)), [RGB(0,0,1), RGB(0,0,0)], lchoices=range(20, stop=70, length=15))
for (i, worm) in enumerate(all_worms)
    lines!(ax_3d, worm.positions_cm, color=colors_traj[i], linewidth=4)
end
display(GLMakie.Screen(), fig_3d_traj) 

interactive_fig = view_binary_stack_with_tracks_glmakie(binary_stack, all_worms[end] )

