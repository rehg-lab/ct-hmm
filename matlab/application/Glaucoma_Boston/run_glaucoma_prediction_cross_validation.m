%clear all; close all; clc

%% setup output directory
top_out_dir = './output_glaucoma_models_cv10';
mkdir(top_out_dir);
global fp_log;
str = sprintf('%s/log.txt', top_out_dir);
fp_log = fopen(str, 'wt');

%% load eye data from the datasheet
num_min_hist_visit = 4; % 4 history visits
min_num_visit = num_min_hist_visit + 1; % at least one future visit
is_parse_glaucoma_data = 0;

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
    temp = load('./data_orig/Boston_qual_eye_list');
    eye_seq_list = temp.qual_eye_list;
end 

fclose(fp_log);

%% data setting
global data_setting;
data_setting = struct;
%data_setting = 0; % for reset
data_setting.dim = 2;
data_setting.type_name_ls = {'VFI', 'RNFL'};
data_setting.dim_value_range_ls{1} = [100.5 99.5 98 96 93 90 85 80:(-10):20]; % dim1
data_setting.dim_value_range_ls{2} = [130:(-5):80 70:(-10):30]; % dim2
data_setting.draw_origin_lefttop = 1;

%% cross validation setting
% num_cv = 10;
% cv_idx_list = [1:1:10];

%num_cv = 1;
%cv_idx_list = [1:1:1];

num_cv = 2;
cv_idx_list = [1:1:2];

%% start training CT-HMMs: 10 fold cross validation
learn_method = 3; % 1: expm 2: unif 3: eigen
is_outer_soft = 1; % 1: soft, 2: hard
func_glaucoma_train_cross_validation(top_out_dir, eye_seq_list, learn_method, is_outer_soft, num_cv, cv_idx_list);

%% prediction test using the pretrained CT-HMM models
pred_out_dir = top_out_dir;
pretrain_CTHMM_dir = top_out_dir;
func_glaucoma_predict_cross_validation(eye_seq_list, pred_out_dir, pretrain_CTHMM_dir, num_cv);
