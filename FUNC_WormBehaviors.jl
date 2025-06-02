module FUNC_WormBehaviors

using LinearAlgebra, GeometryBasics, Statistics
export compute_behaviors

"""
Compute worm behaviors based on trajectory, shape, and arena geometry.
Returns a vector of symbols (:turning, :pausing, :toward, :away, :along) for each point.
"""
function compute_behaviors(positions, times_s,
                            major_axes, minor_axes,
                            img_center, circle_radius)
    n = length(positions)
    behaviors = Vector{Symbol}(undef, n)
    speeds = zeros(Float32, n)
    # compute speeds between consecutive points
    for i in 2:n
        p1, p2 = positions[i-1], positions[i]
        dt = times_s[i] - times_s[i-1]
        speeds[i] = dt > 0 ? norm(Point2f(p2[1]-p1[1], p2[2]-p1[2])) / dt : 0
    end
    mean_speed = mean(speeds[2:end])

    for i in 1:n
        # detect turning via angle change and circularity
        turning = false
        if i > 1 && i < n
            v1 = positions[i] - positions[i-1]
            v2 = positions[i+1] - positions[i]
            ang = acos(clamp(dot(v1, v2) / ((norm(v1)*norm(v2) + eps())), -1, 1))
            circ_ratio = minor_axes[i] / major_axes[i]
            turning = ang > π/8 && circ_ratio > 0.8
        end
        # detect pausing by low speed
        pausing = speeds[i] < 0.1f0 * mean_speed
        # radial velocity
        rv = 0f0
        if i < n
            r1 = Point2f(positions[i][1], positions[i][2]) - img_center
            r2 = Point2f(positions[i+1][1], positions[i+1][2]) - img_center
            rv = norm(r2) - norm(r1)
        end
        # assign behavior
        if turning
            behaviors[i] = :turning
        elseif pausing
            behaviors[i] = :pausing
        elseif rv > 0.01
            behaviors[i] = :away
        elseif rv < 0.01
            behaviors[i] = :toward
        else
            behaviors[i] = :along
        end
    end
    return behaviors
end

end # module
