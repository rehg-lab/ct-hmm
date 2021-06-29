close all; clc

%% please put your datasheet in datasheet_folder 
datasheet_folder = 'Glaucoma_Boston_parse';
addpath(datasheet_folder);

%% specifiy the datasheet name here
%datasheet_name = 'yuying2_cal.xlsx'; % this is the datasheet that Igor sent to me in Aug, 2016 (RNFL already calibrated)
datasheet_name = 'Pitt_cthmm_116_6m.xlsx';
datasheet_file = sprintf('%s/%s', datasheet_folder, datasheet_name);

%% load eye data from the datasheet (here, we include all eyes with at least 5 visit data)
global num_min_hist_visit;
num_min_hist_visit = 4; % 4 history visits
min_num_visit = num_min_hist_visit + 1; % at least one future visit, this means we will include all eyes which have at least 5 visits
is_parse_glaucoma_data = 1;

%% setup output directory
top_out_dir = sprintf('output_predict_cv10_minvisit%d_PittNew', min_num_visit);
mkdir(top_out_dir);

global fp_log;
str = sprintf('%s/log.txt', top_out_dir);
fp_log = fopen(str, 'wt');

data_parse_output_dir = sprintf('%s/data_stat', top_out_dir);
mkdir(data_parse_output_dir);

%% parse glaucoma Boston dataset
if (is_parse_glaucoma_data == 1)
    %% select patients
    is_extract_g_gs = 1; % extract glaucoma and glaucoma suspect      
    disp('Start parsing glaucoma datasheet...');
    
    [all_eye_list, qual_eye_list, num_visit_list, trace_year_list, data1_list, data2_list, baseline_age_list, age_list] = ...
        func_parse_glaucoma_Boston_datasheet_calibrated_general(datasheet_file, data_parse_output_dir, is_extract_g_gs, min_num_visit);
		
    %% select training group
    age_l = 0; age_u = 120; % age low bound and upper bound, select all ages
    [eye_seq_list, select_trace_year_list, select_num_visit_list] = func_glaucoma_select_train_visit_age(age_l, age_u, qual_eye_list);        
    save('qual_eye_list', 'qual_eye_list');
else
    temp = load('qual_eye_list');
    eye_seq_list = temp.qual_eye_list;
end 

%% data setting
global data_setting;
data_setting = struct;
data_setting.dim = 2;
data_setting.type_name_ls = {'VFI', 'RNFL'};
%% for visualization, draw data origin at the left top corner
data_setting.draw_origin_lefttop = 1;

%% setup state grid
data_setting.dim_value_range_ls{1} = [100 99 98 96 94 92 90 85 80 75 70:(-10):0]; % dim1
data_setting.dim_value_range_ls{2} = [110:(-5):30]; % dim2

%% start training CT-HMMs: 
% 10 fold cross validation
num_cv = 2;
cv_idx_list = [1:1:2];

% num_cv = 10;
% cv_idx_list = [1:1:10];

learn_method = 4; % 1: expm, 2: unif, 3: eigen, 4: hybrid (expm in the first iteration, then eigen for the rest)
is_outer_soft = 1; % 1: soft, 2: hard
% func_glaucoma_train_cross_validation(top_out_dir, eye_seq_list, learn_method, is_outer_soft, num_cv, cv_idx_list);

%% prediction test using the pre-trained CT-HMM models
pred_out_dir = top_out_dir;
pretrain_CTHMM_dir = top_out_dir;

func_glaucoma_predict_cross_validation(eye_seq_list, pred_out_dir, pretrain_CTHMM_dir, num_cv);
