% setup the folder for the output figures (./ means the script's current folder)
output_dir = './';  

% the learned model you want to visualize
model_dir = 'application\Glaucoma_Boston\output_vis_minvisit5_PittNew_tol10-6_def7_0423\age_0_120\CV_1\Iter_64';

model_path = sprintf('%s/learn_variables', model_dir);
load(model_path);

addpath('learn');
addpath('MD\MD_vis\2D_vis');
addpath('MD\MD_vis');


%% state dwelling drawing color
global dwelling_time_draw_range_list;
global dwelling_time_draw_color_list;

%dwelling_time_draw_range_list = [0 2 5 100];
%dwelling_time_draw_color_list = [255 0 0 ; 255 0 0; 0 255 0; 0 255 0];

dwelling_time_draw_range_list = [0 10 100];
%dwelling_time_draw_color_list = [255 255 255; 0 0 0; 0 0 0]; % white to
%black
dwelling_time_draw_color_list = [0 0 0; 255 255 255; 255 255 255]; % black to white

CTHMM_learn_vis_Q_mat(output_dir);

