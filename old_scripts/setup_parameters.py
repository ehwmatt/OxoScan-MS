def setup_params(minimum_retention_time=0,maximum_retention_time=15,retention_time_step=1.1/60,retention_time_merging_distance=8,minimum_precursor_mz=800,maximum_precursor_mz=1600,precursor_mz_step=2,precursor_merging_distance=2,grid_unit_size=8,npeaks_for_persistance=1000,selected_diagnostic_ion='138.055',data_folder='/camp/lab/ralserm/working/Glycoproteomics/2020-Project3-dev-python/DataFromDIANN'):
    
    global min_rt
    min_rt = minimum_retention_time
    global max_rt
    max_rt = maximum_retention_time
    global step_rt
    step_rt = retention_time_step
    global merge_distance_rt
    merge_distance_rt = retention_time_merging_distance
    
    global min_mz
    min_mz = minimum_precursor_mz
    global max_mz
    max_mz = maximum_precursor_mz
    global step_mz
    step_mz = precursor_mz_step
    global merge_distance_mz
    merge_distance_mz = precursor_merging_distance
    
    global grid_units
    grid_units = grid_unit_size
    global npeaks
    npeaks = npeaks_for_persistance

    global selected_ion
    selected_ion = selected_diagnostic_ion
    
    global file_extension_ms
    file_extension_ms = ".wiff.dia.extracted.txt"
    global file_extension_peaks
    file_extension_peaks = ".peaks.csv"
    global file_extension_simple_peaks
    file_extension_simple_peaks = ".simple_peaks.csv"

    global central_peaks_filename
    central_peaks_filename = data_folder + "/central_peaks.csv"
    global metatable_filename
    metatable_filename = data_folder + "/metatable.tsv"
    global features_filename
    features_filename = data_folder + "/features.csv"
    global simple_features_filename
    simple_features_filename = data_folder + "/features_simple.csv"
    
    print('Parameters for analysis have been updated')