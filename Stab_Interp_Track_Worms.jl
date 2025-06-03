# filepath: Stab_Interp_Track_Worms.jl
module Stab_Interp_Track_Worms

using TiffImages, GLMakie, Observables
using LinearAlgebra, Statistics, Colors, GeometryBasics
using ProgressMeter, Base.Threads
using ImageMorphology

include("FUNC_ImageProcessing.jl"); using .FUNC_ImageProcessing
include("FUNC_WormFinder.jl"); using .FUNC_WormFinder
include("FUNC_ObjectFind.jl"); using .FUNC_ObjectFind
include("FUNC_WormBehaviors.jl"); using .FUNC_WormBehaviors

include("Plotting_Aux.jl"); using .Plotting_Aux: plot_lightdark

# Conversion constants
const CM_PER_PIXEL = 1.0 / 58.6
const MSEC_PER_FRAME = 50.0 * 2.0
const MIN_AREA = 15
const MAX_AREA = 500
const MAJOR_MINOR_MAX_RATIO = 10.0
const WORM_DISTANCE_SEARCHADD = 30.0

"""
Process a list of TIFF files, track worm positions, and optionally plot 2D light/dark trajectories.
# Arguments
- filepaths::Vector{String}: paths to TIFF stacks
- thresholds::Vector{Float64}: segmentation thresholds (length 1 or matches filepaths)
- do_plot_lightdark::Bool=false: whether to call `plot_lightdark` on results
# Returns
- results::Vector{WormData}
- bitstacks::Vector{BitArray{3}}
"""
function track_worms(filepaths::Vector{String}; thresholds::Vector{Float64}=fill(1.05, length(filepaths)), do_plot_lightdark::Bool=false)
    n = length(filepaths)
    if length(thresholds) == 1
        thresholds = fill(thresholds[1], n)
    elseif length(thresholds) != n
        error("`thresholds` length must be 1 or equal to number of filepaths ($n); got $(length(thresholds)).")
    end
    results = Vector{WormData}(undef, n)
    bitstacks = Vector{BitArray{3}}(undef, n)
    
    @threads for idx in 1:n
        filepath = filepaths[idx]
        thr = thresholds[idx]
        println("Processing: $filepath with threshold $thr on thread $(threadid())")
        raw = TiffImages.load(filepath)
        img_rows, img_cols, _ = size(raw)
        data = 1.0 .- Float32.(raw)
        ref = median(data, dims=3)
        ds = data ./ (ref .+ eps(Float32(1.0)))
        mask = dilate(erode(ds .> thr))
        local_stack = remove_outlier_objects(mask, MIN_AREA, MAX_AREA, MAJOR_MINOR_MAX_RATIO)
        pos, areas, maj, minr = create_track(local_stack, true, WORM_DISTANCE_SEARCHADD)
        beh = compute_behaviors(pos, [p[3]*MSEC_PER_FRAME/1000 for p in pos], maj, minr, Point2f(img_cols/2, img_rows/2), min(img_rows, img_cols)/2)
        times_s = [p[3]*MSEC_PER_FRAME/1000 for p in pos]
        pos_cm = [Point3f(p[1]*CM_PER_PIXEL, (img_rows-p[2])*CM_PER_PIXEL, times_s[i]) for (i,p) in enumerate(pos)]
        results[idx] = WormData(filepath, pos, areas, maj, minr, img_rows, img_cols, times_s, pos_cm, beh)
        bitstacks[idx] = local_stack
    end

    if do_plot_lightdark
        plot_lightdark(results)
    end

    return results, bitstacks
end

export track_worms

end # module
