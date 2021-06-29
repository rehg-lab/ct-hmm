%clear all; 
close all;
clc

addpath('../');
addpath('../../prediction');
addpath('../../decode');
addpath('../../decode/baseline_SSA');
addpath('../../learn');
addpath('../../common');
addpath('../../MD/MD_model');
addpath('../../MD/MD_vis');
addpath('../../MD/MD_vis/2D_vis');
addpath('Glaucoma_Boston_parse');

CTHMM_model_dir = 'C:/exp_output/output_glaucoma_vis_minvisit5_hybrid_statedef2/age_0_110/CV_1/Iter_18';
CTHMM_learn_load_para(CTHMM_model_dir);

global is_draw_learn_Q_mat;
is_draw_learn_Q_mat = 1;

%% state dwelling drawing color
global dwelling_time_draw_range_list;
global dwelling_time_draw_color_list;
dwelling_time_draw_range_list = [0 2 5 100];
dwelling_time_draw_color_list = [255 0 0 ; 255 0 0; 0 255 0; 0 255 0];

top_out_folder = 'output_vis';

%% draw learning figures
%CTHMM_learn_vis_Q_mat(top_out_folder);
%is_vis_nij = 1; % vis nij
%vis_threshold = 0.1;
%is_draw_text = 1.0;
%CTHMM_vis_2D_Q_mat(top_out_folder, is_vis_nij, vis_threshold, is_draw_text);

%vis_mat_type = 1; % nij
%vis_threshold = 0.000001;
%vis_threshold = 0.0;
%is_draw_text = 1;
%[saved_fig] = CTHMM_vis_2D_Q_mat_new(top_out_folder, vis_mat_type, vis_threshold, is_draw_text);

%vis_mat_type = 3; % vij
%vis_threshold = 0.000001;
%is_draw_text = 1;
%[saved_fig] = CTHMM_vis_2D_Q_mat_new(top_out_folder, vis_mat_type, vis_threshold, is_draw_text);


%% draw future prediction 
%% 2, 5, 10, 20 years later

%out_dir = sprintf('%s/future_state', top_out_folder);
%mkdir(out_dir);
%vis_threshold = 0.000001;
%for future_year = 5:5:50
%    future_month = future_year * 12.0;
%    CTHMM_vis_2D_most_probable_future_state(out_dir, vis_threshold, future_month);
%end

out_dir = sprintf('%s/future_path_20years', top_out_folder);
mkdir(out_dir);
state_count_threshold = 10;
future_month = 20 * 12.0; % future 50 years
CTHMM_vis_2D_most_probable_future_path(out_dir, state_count_threshold, future_month);
