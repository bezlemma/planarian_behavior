module Stab_Interp_Track_Worms_LightDark

# Dependencies
using TiffImages, GLMakie, Observables
using LinearAlgebra, Statistics, Colors, GeometryBasics
using ProgressMeter, Base.Threads
using ImageMorphology

# Internal utilities
include("FUNC_ImageProcessing.jl"); using .FUNC_ImageProcessing
include("FUNC_WormFinder.jl"); using .FUNC_WormFinder
include("FUNC_ObjectFind.jl"); using .FUNC_ObjectFind
include("FUNC_WormBehaviors.jl"); using .FUNC_WormBehaviors

# Constants
const CM_PER_PIXEL = 1.0 / 58.6
const MSEC_PER_FRAME = 50.0 * 2.0
const MIN_AREA = 15
const MAX_AREA = 500
const MAJOR_MINOR_MAX_RATIO = 10.0
const WORM_DISTANCE_SEARCHADD = 30.0

# Main tracking function
function track_worms(filepaths::Vector{String}; thresholds::Vector{Float64}=fill(1.05, length(filepaths)), do_plot_lightdark::Bool=true)
    n = length(filepaths)

    # Normalize thresholds: allow single value or matching number of files
    if length(thresholds) == 1
        thresholds = fill(thresholds[1], n)
    elseif length(thresholds) != n
        error("`thresholds` length must be 1 or equal to number of filepaths ($n); got $(length(thresholds)).")
    end
    results = Vector{WormData}(undef, n)
    bitstacks = Vector{BitArray{3}}(undef, n)

    @threads for idx in 1:n
        filepath = filepaths[idx]
        local_threshold = thresholds[idx]
        println("Processing: $filepath with threshold $local_threshold on thread $(threadid())")
        raw_data = TiffImages.load(filepath)
        img_rows, img_cols, num_frames = size(raw_data)

        data = 1.0 .- Float32.(raw_data)
        reference = median(data, dims=3)
        divided_stack = data ./ (reference .+ eps(Float32(1.0)))
        binary_mask = dilate(erode(divided_stack .> local_threshold))
        local_stack = remove_outlier_objects(binary_mask, MIN_AREA, MAX_AREA, MAJOR_MINOR_MAX_RATIO)

        positions, areas_2, major_axes, minor_axes = create_track(local_stack, true, WORM_DISTANCE_SEARCHADD)
        behaviors = compute_behaviors(positions,
                                     [p[3] * MSEC_PER_FRAME / 1000.0 for p in positions],
                                     major_axes, minor_axes,
                                     Point2f(img_cols/2, img_rows/2),
                                     min(img_rows, img_cols)/2)

        times_s = [p[3] * MSEC_PER_FRAME / 1000.0 for p in positions]
        pos_cm = [Point3f(p[1]*CM_PER_PIXEL, (img_rows - p[2])*CM_PER_PIXEL, times_s[i]) for (i,p) in enumerate(positions)]

        results[idx] = WormData(filepath, positions, areas_2, major_axes, minor_axes, img_rows, img_cols, times_s, pos_cm, behaviors)
        bitstacks[idx] = local_stack
    end

    return results, bitstacks
end

# Export API
export track_worms

end # module
