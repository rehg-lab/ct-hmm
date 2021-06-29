%% decode the glaucoma dataset
clear all; close all; clc

%% add decoding paths
addpath('../../decode/baseline_SSA');
addpath('../../decode');

%% setup output directory
top_out_dir = '../../../output_glaucoma_decoding';
mkdir(top_out_dir);

global fp_log;
str = sprintf('%s/log.txt', top_out_dir);
fp_log = fopen(str, 'wt');

%% load eye data from the datasheet
num_min_hist_visit = 4; % 4 history visits
min_num_visit = num_min_hist_visit + 1; % at least one future visit
is_parse_glaucoma_data = 1;

%% parse glaucoma Boston dataset
if (is_parse_glaucoma_data == 1)
    %% select patients
    is_extract_g_gs = 1; % extract glaucoma and glaucoma suspect  
    str_Dx = 'g/gs'; % diagnosis code
    dataset_name = 'Boston';
    datasheet_folder = 'Glaucoma_Boston_parse';       
    addpath(datasheet_folder);    
    
    [all_eye_list, qual_eye_list, num_visit_list, trace_year_list, VFI_list, RNFL_list, baseline_age_list, age_list] = func_parse_glaucoma_Boston_datasheet(datasheet_folder, is_extract_g_gs, min_num_visit);
    
        
    %% select training group
    age_l = 0; age_u = 120; % age low bound and upper bound, select all ages
    [eye_seq_list] = func_glaucoma_select_train_visit_age(age_l, age_u, qual_eye_list);        
    save('Boston_qual_eye_list', 'qual_eye_list');    
else
    temp = load('Boston_qual_eye_list');
    eye_seq_list = temp.qual_eye_list;
end 

%% data setting
global data_setting;


data_setting.dim = 2;
data_setting.type_name_ls = {'VFI', 'RNFL'};
%data_setting.dim_value_range_ls{1} = [100.5 99.5 98 96 93 90 85 80:(-10):20]; % dim1
%data_setting.dim_value_range_ls{2} = [130:(-5):80 70:(-10):30]; % dim2
data_setting.dim_value_range_ls{1} = [100.5 99.5 98 96 93 90 85 80:(-10):20]; % dim1: VFI
data_setting.dim_value_range_ls{2} = [130:(-5):50 40 30]; % dim2: RNFL
data_setting.draw_origin_lefttop = 1;

%% cross validation setting
num_cv = 1;
cv_idx_list = [1];

%% start training CT-HMM
global is_draw_learn_Q_mat;
is_draw_learn_Q_mat = 1;

learn_method = 3; % eigen
is_outer_soft = 1; % soft
func_glaucoma_train_cross_validation(top_out_dir, eye_seq_list, learn_method, is_outer_soft, num_cv, cv_idx_list);
                        
%% load the trained CT-HMM model
cv_idx = 1;
str = sprintf('%s/CV_%d/num_iter.txt', top_out_dir, cv_idx);    
%fp = fopen(str, 'rt');
%num_iter = fscanf(fp, '%d');
num_iter = 20;

CTHMM_model_dir = sprintf('%s/CV_%d/Iter_%d', top_out_dir, cv_idx, num_iter);    
CTHMM_learn_load_para(CTHMM_model_dir);  % now the Q_mat is loaded

%% start decoding every sequence
num_obs_seq = size(eye_seq_list, 1);
global is_decode_simulation;
is_decode_simulation = 0;
global is_glaucoma_decoding;
is_glaucoma_decoding = 1;

for i = 1:num_obs_seq    
    
    i
    
    %% continuous-time decoding
    eye_seq_list{i}.has_compute_data_emiss_prob = 0;    
    tic;
    expect_dur_method = 1; % expm
    [decoded_conti_state_seq, decoded_conti_dwelltime_seq, decoded_obs_state_seq] = CTHMM_decode_complete_path_viterbiexpect(eye_seq_list{i}, expect_dur_method);
    tEnd = toc;
    
    
    %decoded_conti_state_seq = [2 1 60 35]';
    %decoded_conti_dwelltime_seq = [32.2783  6.8847 2.2632 2.5738]';
    %decoded_obs_state_seq = [2 2 2 1 35]';
       
    
    %% visualize the decoded sequence with nij count
    %is_vis_nij = 1;
    %vis_threshold = 0.1;
    %is_draw_text = 0.0;
    %[saved_fig] = CTHMM_vis_2D_Q_mat(top_out_dir, is_vis_nij, vis_threshold, is_draw_text);
    
    %% draw continuous decoded  trajectory
    
    is_load_ave_trend_fig = 0;
    %load_figname = saved_fig;
    load_figname = '';
    out_filename = sprintf('%s/decoded_traj_%03d', top_out_dir, i);
    [saved_fig] = CTHMM_vis_2D_decoded_conti_traj(eye_seq_list{i}, decoded_conti_state_seq, decoded_conti_dwelltime_seq, decoded_obs_state_seq, out_filename, is_load_ave_trend_fig, load_figname);     
            
    %close(gcf);
    
end

%str = sprintf('CTHMM_decode_batch: %d min, %f sec\n', floor(tEnd/60),rem(tEnd,60));
%CTHMM_print_log(str);                               

fclose(fp_log);
