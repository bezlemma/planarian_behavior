module FUNC_ImageProcessing

using AbstractFFTs
using FFTW
using Base.Threads
using ProgressMeter
using Statistics  

export create_border_image, fftshift2d, phase_correlate, translate_image, stabilize

function create_border_image(img::AbstractMatrix{T}) where T
    rows, cols = size(img)
    border_only_img = copy(img)
    center_h = floor(Int, rows * 0.9)
    center_w = floor(Int, cols * 0.9)
    if center_h < rows && center_w < cols && center_h > 0 && center_w > 0
        margin_y = rows - center_h; offset_y = margin_y ÷ 2; start_y_center = offset_y + 1
        end_y_center = start_y_center + center_h - 1
        margin_x = cols - center_w; offset_x = margin_x ÷ 2; start_x_center = offset_x + 1
        end_x_center = start_x_center + center_w - 1
        actual_start_y = max(1, start_y_center); actual_end_y = min(rows, end_y_center)
        actual_start_x = max(1, start_x_center); actual_end_x = min(cols, end_x_center)
        if actual_start_y <= actual_end_y && actual_start_x <= actual_end_x
            border_only_img[actual_start_y:actual_end_y, actual_start_x:actual_end_x] .= T(0)
        end
    elseif center_h >= rows || center_w >= cols
        border_only_img .= T(0)
    end
    return border_only_img
end

function fftshift2d(mat)
    sy, sx = size(mat)
    return circshift(mat, (sy ÷ 2, sx ÷ 2))
end

function phase_correlate(ref_img_orig, moving_img_orig)
    ref_img_cropped, moving_img_cropped = if size(ref_img_orig) != size(moving_img_orig)
            r_h, r_w = size(ref_img_orig); m_h, m_w = size(moving_img_orig)
            min_h = min(r_h, m_h); min_w = min(r_w, m_w)
            ref_img_orig[1:min_h, 1:min_w], moving_img_orig[1:min_h, 1:min_w]
        else
            ref_img_orig, moving_img_orig
        end
    ref_img_border = create_border_image(ref_img_cropped)
    moving_img_border = create_border_image(moving_img_cropped)
    fft_ref = fft(Float32.(ref_img_border)); fft_moving = fft(Float32.(moving_img_border))
    cross_power_spectrum = fft_ref .* conj(fft_moving)
    magnitude = abs.(cross_power_spectrum)
    epsilon = Float32(1e-10)
    normalized_cross_power_spectrum = cross_power_spectrum ./ (magnitude .+ epsilon)
    correlation_matrix = fftshift2d(real.(ifft(normalized_cross_power_spectrum)))
    _max_val, max_idx_linear = findmax(correlation_matrix)
    max_idx_cartesian = CartesianIndices(size(correlation_matrix))[max_idx_linear]
    center_y, center_x = (size(correlation_matrix) .÷ 2) .+ 1
    dy = max_idx_cartesian[1] - center_y; dx = max_idx_cartesian[2] - center_x
    return dy, dx
end

function translate_image(img::AbstractMatrix{T}, dy, dx, fill_value::T) where T
    rows, cols = size(img)
    translated_img = fill(fill_value, rows, cols)
    src_y_start = max(1, 1 - dy); src_y_end = min(rows, rows - dy)
    src_x_start = max(1, 1 - dx); src_x_end = min(cols, cols - dx)
    dest_y_start = max(1, 1 + dy); dest_y_end = min(rows, rows + dy)
    dest_x_start = max(1, 1 + dx); dest_x_end = min(cols, cols + dx)
    if (src_y_start <= src_y_end) && (src_x_start <= src_x_end) && (dest_y_start <= dest_y_end) && (dest_x_start <= dest_x_end)
        src_height = src_y_end - src_y_start + 1; src_width = src_x_end - src_x_start + 1
        dest_height = dest_y_end - dest_y_start + 1; dest_width = dest_x_end - dest_x_start + 1
        copy_height = min(src_height, dest_height); copy_width = min(src_width, dest_width)
        src_y_end_clipped = src_y_start + copy_height - 1; src_x_end_clipped = src_x_start + copy_width - 1
        dest_y_end_clipped = dest_y_start + copy_height - 1; dest_x_end_clipped = dest_x_start + copy_width - 1
        if copy_height > 0 && copy_width > 0
            translated_img[dest_y_start:dest_y_end_clipped, dest_x_start:dest_x_end_clipped] = img[src_y_start:src_y_end_clipped, src_x_start:src_x_end_clipped]
        end
    end
    return translated_img
end

function stabilize(data,reference)
    @showprogress Threads.@threads for i in 1:lastindex(data, 3)
        current_frame = data[:, :, i]
        dy, dx = phase_correlate(reference, current_frame)
        data[:, :, i] = translate_image(current_frame, -dy, -dx, 0.0f0)
    end
    return data, median(data, dims=3)
end

end # module