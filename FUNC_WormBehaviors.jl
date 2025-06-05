module FUNC_WormBehaviors

using LinearAlgebra, GeometryBasics, Statistics
export compute_behaviors, smooth_behaviors

"""
Compute worm behaviors based on trajectory, shape, and arena geometry.
Returns a vector of symbols (:turning, :pausing, :toward, :away, :along) for each point.
"""
function compute_behaviors(positions, times_s,
                            major_axes, minor_axes,
                            img_center, circle_radius)
    n = length(positions)
    raw_behaviors = Vector{Symbol}(undef, n)
    speeds = zeros(Float32, n)
    
    # compute speeds between consecutive points
    for i in 2:n
        p1, p2 = positions[i-1], positions[i]
        dt = times_s[i] - times_s[i-1]
        speeds[i] = dt > 0 ? norm(Point2f(p2[1]-p1[1], p2[2]-p1[2])) / dt : 0
    end
    mean_speed = mean(speeds[2:end])

    # first pass: detect raw behaviors
    for i in 1:n
        # detect turning via angle change and circularity
        turning = false
        if i > 1 && i < n
            v1 = positions[i] - positions[i-1]
            v2 = positions[i+1] - positions[i]
            ang = acos(clamp(dot(v1, v2) / ((norm(v1)*norm(v2) + eps())), -1, 1))
            circ_ratio = minor_axes[i] / major_axes[i]
            #TODO: Mka this an argument
            turning = ang > π/16 && circ_ratio > 0.6 #was π/8 and 0.8
        end
        
        # detect pausing by low speed - only truly stationary worms
        pausing = speeds[i] < 0.0005f0  # Absolute minimal threshold
        
        # compute radial and tangential components of movement
        r_vel = 0f0  # radial velocity
        t_vel = 0f0  # tangential velocity
        if i < n
            pos_current = Point2f(positions[i][1], positions[i][2])
            pos_next = Point2f(positions[i+1][1], positions[i+1][2])
            
            r1 = pos_current - img_center
            r2 = pos_next - img_center
            
            # radial velocity (change in distance from center)
            r_vel = norm(r2) - norm(r1)
            
            # tangential velocity (component perpendicular to radius)
            if norm(r1) > 1e-6  # avoid division by zero
                movement_vec = pos_next - pos_current
                radial_unit = normalize(r1)
                tangential_component = movement_vec - dot(movement_vec, radial_unit) * radial_unit
                t_vel = norm(tangential_component)
            end
        end
        
        if turning
            raw_behaviors[i] = :turning
        elseif pausing  # Only classify as pausing if truly stationary
            raw_behaviors[i] = :pausing
        else
            pos_current = Point2f(positions[i][1], positions[i][2])
            if norm(pos_current - img_center) < circle_radius*0.8 # worm inside inner 90% of arena
                    raw_behaviors[i] = r_vel > 0 ? :away : :toward
            else
                if abs(r_vel) > t_vel
                    raw_behaviors[i] = r_vel > 0 ? :away : :toward
                else
                    raw_behaviors[i] = :along
                end
            end
        end
    end
    
    # second pass: smooth behaviors to reduce flickering
    smoothed_behaviors = smooth_behaviors(raw_behaviors)
    
    return smoothed_behaviors
end



"""
Smooth behavior sequence to reduce flickering by applying majority filtering
and filling isolated behavior spikes.
"""
function smooth_behaviors(behaviors::Vector{Symbol})
    n = length(behaviors)
    if n <= 2
        return behaviors
    end
    
    smoothed = copy(behaviors)
    window_size = 20  # window for smoothing
    
    # first pass: majority filter
    for i in 2:(n-1)
        window_start = max(1, i - window_size ÷ 2)
        window_end = min(n, i + window_size ÷ 2)
        
        # count behavior frequencies in window
        behavior_counts = Dict{Symbol, Int}()
        for j in window_start:window_end
            behavior_counts[behaviors[j]] = get(behavior_counts, behaviors[j], 0) + 1
        end
        
        # find most common behavior in window
        max_count = 0
        most_common = behaviors[i]
        for (behavior, count) in behavior_counts
            if count > max_count
                max_count = count
                most_common = behavior
            end
        end
        
        # only change if there's a clear majority (more than half the window)
        if max_count > (window_end - window_start + 1) ÷ 2
            smoothed[i] = most_common
        end
    end
    
    # second pass: fill isolated single-frame behaviors
    for i in 2:(n-1)
        if smoothed[i] != smoothed[i-1] && smoothed[i] != smoothed[i+1] && smoothed[i-1] == smoothed[i+1]
            # isolated behavior: replace with neighboring behavior
            smoothed[i] = smoothed[i-1]
        end
    end
    
    # third pass: handle short sequences (2-3 frames) of behaviors
    i = 1
    while i <= n-2
        if smoothed[i] != smoothed[i+1] && smoothed[i+1] != smoothed[i+2]
            # look for pattern like A-B-A or A-B-C where B is isolated
            if i < n-2 && smoothed[i] == smoothed[i+2]
                smoothed[i+1] = smoothed[i]  # fill A-B-A → A-A-A
                i += 3
            elseif i < n-3 && smoothed[i] == smoothed[i+3]
                # pattern A-B-C-A, fill middle with A
                smoothed[i+1] = smoothed[i]
                smoothed[i+2] = smoothed[i]
                i += 4
            else
                i += 1
            end
        else
            i += 1
        end
    end
    
    return smoothed
end

end # module
