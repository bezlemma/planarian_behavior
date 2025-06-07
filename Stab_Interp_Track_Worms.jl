# filepath: Stab_Interp_Track_Worms.jl
module Stab_Interp_Track_Worms

using TiffImages, GLMakie, Observables
using LinearAlgebra, Statistics, Colors, GeometryBasics
using ProgressMeter, Base.Threads
using ImageMorphology
using StatsBase: percentile

include("FUNC_ImageProcessing.jl"); using .FUNC_ImageProcessing
include("FUNC_WormFinder.jl"); using .FUNC_WormFinder
include("FUNC_ObjectFind.jl"); using .FUNC_ObjectFind
include("FUNC_WormBehaviors.jl"); using .FUNC_WormBehaviors

include("Plotting_Aux.jl"); using .Plotting_Aux: plot_lightdark

# Conversion constants
const CM_PER_PIXEL = 1.0 / 58.6
const MSEC_PER_FRAME = 50.0
const MIN_AREA = 15
const MAX_AREA = 500
const MAJOR_MINOR_MAX_RATIO = 10.0
const WORM_DISTANCE_SEARCHADD = 30.0

"""
Process a list of TIFF files, track worm positions.
# Arguments
- filepaths::Vector{String}: paths to TIFF stacks
- thresholds::Vector{Float64}: segmentation thresholds (length 1 or matches filepaths)
# Returns
- results::Vector{WormData}
- bitstacks::Vector{BitArray{3}}
"""
function track_worms(filepaths::Vector{String}; 
                     pixel_percentages::Vector{Float64}=fill(0.25, length(filepaths))) # Default to 0.25%
    n = length(filepaths)
    local actual_pixel_percentages::Vector{Float64} # Ensure type stability
    if length(pixel_percentages) == 1
        actual_pixel_percentages = fill(pixel_percentages[1], n)
    elseif length(pixel_percentages) != n
        error("`pixel_percentages` length must be 1 or equal to number of filepaths ($n); got $(length(pixel_percentages)).")
    else
        actual_pixel_percentages = pixel_percentages
    end
    results = Vector{WormData}(undef, n)
    bitstacks = Vector{BitArray{3}}(undef, n)
    
    @threads for idx in 1:n
        filepath = filepaths[idx]
        current_pixel_percentage = actual_pixel_percentages[idx]
        println("Processing: $filepath with pixel_percentage $current_pixel_percentage% on thread $(threadid())")
        
        raw_stack = TiffImages.load(filepath)
        img_rows, img_cols, num_frames = size(raw_stack)
        
        raw_median_ref_2D = Float32.(median(Float32.(raw_stack), dims=3)[:,:,1])

        local_stack_frames = BitArray{3}(undef, img_rows, img_cols, num_frames)

        for frame_idx in 1:num_frames
            current_raw_frame = Float32.(raw_stack[:,:,frame_idx])
            
            normalized_frame = current_raw_frame ./ (raw_median_ref_2D .+ eps(Float32))
        
            inverted_normalized_frame = 1.0f0 .- normalized_frame
            
            processed_frame = clamp.(inverted_normalized_frame, 0.0f0, 1.0f0)
            
            threshold_value = percentile(vec(processed_frame), 100.0 - current_pixel_percentage)
            
            binary_frame = processed_frame .> threshold_value
            local_stack_frames[:,:,frame_idx] = binary_frame
        end

        local_stack = remove_outlier_objects(local_stack_frames, MIN_AREA, MAX_AREA, MAJOR_MINOR_MAX_RATIO)
        
        pos, areas, maj, minr = create_track(local_stack, true, WORM_DISTANCE_SEARCHADD)
        beh = compute_behaviors(pos, [p[3]*MSEC_PER_FRAME/1000 for p in pos], maj, minr, Point2f(img_cols/2, img_rows/2), min(img_rows, img_cols)/2)
        times_s = [p[3]*MSEC_PER_FRAME/1000 for p in pos]
        pos_cm = [Point3f(p[1]*CM_PER_PIXEL, (img_rows-p[2])*CM_PER_PIXEL, times_s[i]) for (i,p) in enumerate(pos)]
        results[idx] = WormData(filepath, pos, areas, maj, minr, img_rows, img_cols, times_s, pos_cm, beh)
        bitstacks[idx] = local_stack
    end

    return results, bitstacks
end

export track_worms

end # module
