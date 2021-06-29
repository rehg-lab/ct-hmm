For dealing with new dataset, you will need to:

(1) parse the datasheet, and fill a cell array that contains all subjects' longitudinal data, like this

subject_list{s}.num_visit
subject_list{s}.visit_time_list(v) => v for each visit
subject_list{s}.visit_data_list(v, d) => d for each dimension
and the following two repetitive structures for backward compatible for my old visualization related codes:
subject_list{s}.visit_data{v}.time 
subject_list{s}.visit_data{v}.data(d) => d for each dimension

(2) for multi-dimensional disease model, set up the global variable "data_setting" for assigning the grid for each dimension, name for each dimension, etc. (see run_Alzheimer.m as an example)

(3) for visualization, set up global variables "dwelling_time_draw_range_list", "dwelling_time_draw_color_list" (see run_Alzheimer.m)

(4) to visualize learned Q mat, you may reference/use "CTHMM_vis_2D_Q_mat", or "CTHMM_vis_3D_Q_mat" function (under MD/MD_vis folder). 