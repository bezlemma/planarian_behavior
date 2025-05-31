#You give this script a folder path and it will move all tiff files in that folder to subfolders based on the  file name.

using Glob

folder_path = "/Users/dl0346/Documents/PlanarianVideos/May28/TrailFollow_AllConditions_DONE"

for filepath in glob("*.tiff", folder_path)
    filename_with_ext = basename(filepath)
    filename_no_ext, ext = splitext(filename_with_ext)

    subfolder_name = ""
    if length(filename_no_ext) > 4 && all(isdigit, filename_no_ext[end-3:end])
        subfolder_name = filename_no_ext[1:end-4]
        if endswith(subfolder_name, '_')
            subfolder_name = rstrip(subfolder_name, '_')
        end

    end


    target_subfolder = joinpath(folder_path, subfolder_name)
    mkpath(target_subfolder) 

    destination_filepath = joinpath(target_subfolder, filename_with_ext)
    mv(filepath, destination_filepath; force=true) 
end

