%% run all
% clear all; close all; clc

%% load datasheet & normalize data
datasheet_folder = 'LongData';

addpath('../../learn');
addpath('../../decode');
addpath('../../precompute');
addpath('../../likelihood');
addpath('../../common');
addpath('../../common/fast_expm_A_t');
addpath('../../MD/MD_model/');
addpath('../../MD/MD_vis/3D_vis');
addpath('../../MD/MD_vis');

%% set up output directory
global out_dir;
out_dir = 'Output_3D_ABeta_LeftHippo_ADAS_expmsoft_1d4';
mkdir(out_dir);
global fp_log;

global learn_is_know_ground_truth_Q_mat; % we know ground truth as this is simulation test
learn_is_know_ground_truth_Q_mat = 0;
global is_use_distinct_time_grouping; % group visits by time to save computation 
is_use_distinct_time_grouping = 1;
global is_use_individual_Q_mat; % set off, as we now didn't model covariates
is_use_individual_Q_mat = 0;
global is_draw_learn_Q_mat; % we draw Q mat for glaucoma and alzhimer's disease to see the state transition trends\
is_draw_learn_Q_mat = 1; 

%% 3D data setting
global data_setting;
dim = 3;
data_setting.dim = dim;
data_setting.draw_origin_lefttop = 0;
data_setting.max_data_count = 5000;
data_setting.type_name_ls = {'A-Beta level', 'Left Hippocampus', 'ADAS-cog'};
data_setting.dim_value_range_ls{1} = [350:-25:0]; % X: ABeta:0-350   
%data_setting.dim_value_range_ls{1} = [0:25:350 500]; % X: tau:0-500
data_setting.dim_value_range_ls{2} = [5500 4500:-250:1000]; % Y: lefthippo
data_setting.dim_value_range_ls{3} = [0:5:70]; % Z: ADAS: 0-70
use_cognition = 1;
use_brainMRI = 1;
use_biochemical = 1;
%data_idx = [2 4 1];
data_idx = [3 4 1];

global train_group_file_str;
global train_group_title;
train_group_file_str = '3D_ABeta_LeftHippo_ADAScog';
train_group_title = 'ABeta_LeftHippo_ADAScog';

global dwelling_time_draw_range_list;
global dwelling_time_draw_color_list;
dwelling_time_draw_range_list = [0 4 10];
dwelling_time_draw_color_list = [255 0 0; 0 255 0; 0 255 0];
global neighbor_link_setting;
neighbor_link_setting = [1 1 1]; % forwarding links at every dimension

%% create log file
str = sprintf('%s/log.txt', out_dir);
fp_log = fopen(str, 'wt');

is_parse_alzheimer_data = 1;
if (is_parse_alzheimer_data == 1)
    [all_subject_list_Alzheimer] = parse_ADNI1_data(datasheet_folder);
    save('all_subject_list_Alzheimer', 'all_subject_list_Alzheimer');
else
    load('all_subject_list_Alzheimer');
end

%% start creating the state structure and Q mat
global obs_seq_list;
[obs_seq_list] = func_alzheimer_select_ND_subject_data(all_subject_list_Alzheimer, use_cognition, use_brainMRI, use_biochemical, data_idx);    
num_subject = size(obs_seq_list, 1);
train_idx_list = [1:1:num_subject]';

%% create state list
global state_list;
global is_add_state_straight_path;
is_add_state_straight_path = 1;
state_sigma = 0.25;
CTHMM_MD_create_state_list(train_idx_list, state_sigma);
CTHMM_vis_3D_state_statistics();
str = sprintf('%s/state_statistics.png', out_dir);
saveas(gcf, str, 'png');

%% create state reachibility matrix
CTHMM_precompute_state_reach_mat();

%% create Q structure
CTHMM_MD_create_Q_mat_struct();

%% set up initial state probability distribution    
num_state = size(state_list, 1);
global state_init_prob_list;
state_init_prob_list = zeros(num_state, 1);
state_init_prob_list(:) = 1.0 / num_state;    

%% init learning Q mat
global Q_mat_init;
init_ave_state_dwell = 2; % 2 years

%is_add_random_perturb = 0;
%perturb_amount = 0;    

is_add_random_perturb = 1; % for eigen method
perturb_amount = 0.05;

Q_mat_init = CTHMM_learn_init_Q_mat(init_ave_state_dwell, is_add_random_perturb, perturb_amount); % assign uniform probability to each link
    
%% do viterbi training using a list of subjects
num_train_id = size(train_idx_list, 1);

%% learning
global learn_converge_tol;
learn_converge_tol = 10^(-5);

%learn_method = 1; % expm
%is_outer_soft = 1; % soft

learn_method = 3; % expm, unif, eigen
is_outer_soft = 1; % soft, hard

max_iter = 1000;
CTHMM_learn_para_EM(learn_method, max_iter, is_outer_soft, Q_mat_init, train_idx_list); 
     
%==============================
%% close log file
fclose(fp_log);

%% load parameter and redo visualization
%load_folder = 'Output_3D_ABeta_LeftHippo_ADAS_expmsoft_1d4\\Iter_20';   
%CTHMM_load_learned_para(load_folder);    
%is_draw_text = 0;    
%is_vis_nij = 1;
%ni_vis_threshold = 0.00;
%CTHMM_vis_3D_Q_mat(out_dir, is_vis_nij, ni_vis_threshold);  
