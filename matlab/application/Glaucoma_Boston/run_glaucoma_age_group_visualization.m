clear; close all; clc

%% load eye data from the datasheet
min_num_visit = 2;
is_parse_glaucoma_data = 0;

%% parse glaucoma Boston datasheet
if (is_parse_glaucoma_data == 1)    
    is_extract_g_gs = 1;   % extract glaucoma and glaucoma suspect  
    str_Dx = 'g/gs'; % diagnosis code
    dataset_name = 'Boston';
    datasheet_folder = 'Glaucoma_Boston_parse';       
    addpath(datasheet_folder);
    [all_eye_list, qual_eye_list, num_visit_list, trace_year_list, VFI_list, RNFL_list, baseline_age_list, age_list] = func_parse_glaucoma_Boston_datasheet(datasheet_folder, is_extract_g_gs, min_num_visit);    
    save('Boston_qual_eye_list', 'qual_eye_list');
else
    temp = load('./data_orig/Boston_qual_eye_list');
    qual_eye_list = temp.qual_eye_list;
end 

%% setup output directory
% top_out_dir = '../../../output_glaucoma_10cv_models_expm';
top_out_dir = '../../../cthmm_outputs/model_orig/';
mkdir(top_out_dir);

global fp_log;
log_filename = sprintf('%s/log.txt', top_out_dir);
fp_log = fopen(log_filename, 'wt');

learn_method = 4; % % 1: expm, 2: unif, 3: eigen, 4: hybrid (expm in the first iteration, then eigen for the rest)
is_outer_soft = 1; % soft

learn_method

%% data setting
global data_setting;
%data_setting = 0; % reset
%% 2015/6/28: test visualization for different age groups
data_setting.dim = 2;
data_setting.type_name_ls = {'VFI', 'RNFL'};

%data_setting.dim_value_range_ls{1} = [100.5 99.5 98 96 93 90 85 80:(-10):20]; % dim1: VFI
%data_setting.dim_value_range_ls{2} = [130:(-5):50 40 30]; % dim2: RNFL

%try a sparser setting
%data_setting.dim_value_range_ls{1} = [100.5 99.5 98 96 93 90 85 80:(-5):60 50:(-10):0]; % dim1: VFI
%data_setting.dim_value_range_ls{2} = [140:(-10):110 105:(-5):60 50:(-10):0]; % dim2: RNFL

% for debug, try a more sparse setting
%data_setting.dim_value_range_ls{1} = [100.5 99.5 98 96 93 90 85 80:(-10):0]; % dim1: VFI
%data_setting.dim_value_range_ls{2} = [140:(-10):110 105:(-5):60 50:(-10):0]; % dim2: RNFL

data_setting.dim_value_range_ls{1} = [100.5 99.5 98 96 93 90 85 80:(-10):0]; % dim1: VFI
data_setting.dim_value_range_ls{2} = [140:(-10):0]; % dim2: RNFL

data_setting.draw_origin_lefttop = 1;

num_cv = 1; % in fact, this means just train 1 model
cv_idx_list = 1;

%% select training group
%% all age
age_l = 0;
age_u = 120;
[select_eye_seq_list] = func_glaucoma_select_train_visit_age(age_l, age_u, qual_eye_list);        
func_glaucoma_train_cross_validation(top_out_dir, select_eye_seq_list, learn_method, is_outer_soft, num_cv, cv_idx_list);

% %% age [0-70]
% age_l = 0;
% age_u = 70;
% [select_eye_seq_list] = func_glaucoma_select_train_visit_age(age_l, age_u, qual_eye_list);        
% func_glaucoma_train_cross_validation(top_out_dir, select_eye_seq_list, learn_method, is_outer_soft, num_cv, cv_idx_list);

% %% age [70-120]
% age_l = 70;
% age_u = 120;
% [select_eye_seq_list] = func_glaucoma_select_train_visit_age(age_l, age_u, qual_eye_list);
% func_glaucoma_train_cross_validation(top_out_dir, select_eye_seq_list, learn_method, is_outer_soft, num_cv, cv_idx_list);
