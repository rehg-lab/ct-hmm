%we want to have a colormap map representation of the 2D CT-HMM output.
%For each of the very detailed states (meaning every 2 microns as the data range of each state for both RNFL and GCIPL thickness), we want to show different color representing the "mean" direction of transitions together with the intensity level (brightness) representing the number of observations in the model.
%Let's say there is a state where 100 subjects passed through. Sixty out of 100 stayed there for a unit of time (say 6 months), 20 moved horizontally (GCIPL progression), 12 moved diagonally, and 8 moved vertically (RNFL progression). This means 20% at 0 degree, 12% at 45 degree, and 8 % at 90 degree. This means 40% moved to 31.5 degree (averaged direction). Since 60% did not move at all, the intensity should be 40%. This is normalized to individual state population.
%Alternatively, we can use the entire population to normalize the intensity level.

% setup the folder for the output figures
top_out_folder = 'out_colormap_densedef4_expm';
mkdir(top_out_folder);

%% read the number of included subjects from the log.txt file
%log_file = 'output_vis_minvisit5_PittNew_tol10-6_def7_0423/log.txt';
log_file = 'output_vis_minvisit5_PittNew_tol10-6_def4_0427_2018_expm_rerun/log.txt';

fp = fopen(log_file, 'rt');
num_total_subject = fscanf(fp, 'num included subject = %d');
fclose(fp);

% the learned model you want to visualize
%model_dir = 'output_vis_minvisit5_PittNew_tol10-6_def7_0423/age_0_120/CV_1/Iter_64';
model_dir = 'output_vis_minvisit5_PittNew_tol10-6_def4_0427_2018_expm_rerun/age_0_120/CV_1/Iter_31';
model_path = sprintf('%s/learn_variables', model_dir);
load(model_path);

future_time_list = [6 12 24 36 48 60 72 84 96 108 120]; % these are future times in months (0.5-10 years)

%future_time_list = [6 24 120];
is_overlay_text = 1;

%% thw following two parameters can be used to filter out some noise, and make the figure more clean
min_count_to_draw_state = 0.0; % minimum subject count threshold to decide whether to draw a state or not. Set 0 if you want to draw all the created states. This is used to filter out the noise 
min_accum_dwell_time_to_draw_absorb_state = 0.5; % minimum accumulated dwelling time (in month) for drawing an absorbing state. Here set 1.0 means at least one month accumulated dwelling time from all population
func_vis_glaucoma_2D_Q_mat_colormap(top_out_folder, future_time_list, num_total_subject, is_overlay_text, min_count_to_draw_state, min_accum_dwell_time_to_draw_absorb_state);


is_redraw_Q_mat = 0;

global is_draw_learn_Q_mat;
global dwelling_time_draw_range_list;
global dwelling_time_draw_color_list;

    
if (is_redraw_Q_mat == 1)

    addpath('../../learn');
    addpath('../../MD/MD_vis/2D_vis');
    addpath('../../MD/MD_vis');

    is_draw_learn_Q_mat = 1;    
    dwelling_time_draw_range_list = [0 10 100];
    %dwelling_time_draw_color_list = [255 255 255; 0 0 0; 0 0 0]; % white to
    %black
    dwelling_time_draw_color_list = [0 0 0; 255 255 255; 255 255 255]; % black to white
    %CTHMM_learn_vis_Q_mat(top_out_folder);

    vis_mat_type = 1; % nij
    %vis_threshold = 0.000001;
    vis_threshold = 0.0;
    is_draw_text = 1;
    state_is_vis_tau_i = 1; % draw total time at each state as state size
    CTHMM_vis_2D_Q_mat_new(top_out_folder, vis_mat_type, vis_threshold, is_draw_text, state_is_vis_tau_i);


    vis_mat_type = 1; % nij
    %vis_threshold = 0.000001;
    vis_threshold = 0.0;
    is_draw_text = 1;
    state_is_vis_tau_i = 0; % draw n_i = \sum n_ij, as state size
    CTHMM_vis_2D_Q_mat_new(top_out_folder, vis_mat_type, vis_threshold, is_draw_text, state_is_vis_tau_i);

    vis_mat_type = 3; % vij
    %vis_threshold = 0.000001;
    vis_threshold = 0.0;
    is_draw_text = 1;
    state_is_vis_tau_i = 0; % draw n_i = \sum n_ij, as state size
    CTHMM_vis_2D_Q_mat_new(top_out_folder, vis_mat_type, vis_threshold, is_draw_text, state_is_vis_tau_i);
 
end