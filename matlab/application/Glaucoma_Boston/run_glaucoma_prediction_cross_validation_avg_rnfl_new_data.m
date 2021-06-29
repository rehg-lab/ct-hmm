close all; clc

%% please put your datasheet in datasheet_folder 
datasheet_folder = 'Glaucoma_Boston_parse';
addpath(datasheet_folder);

%% specifiy the datasheet name here
%datasheet_name = 'yuying2_cal.xlsx'; % this is the datasheet that Igor sent to me in Aug, 2016 (RNFL already calibrated)
% datasheet_name = 'Pitt_cthmm_116_6m.xlsx';
% datasheet_file = sprintf('%s/%s', datasheet_folder, datasheet_name);

%% load eye data from the datasheet (here, we include all eyes with at least 5 visit data)
global num_min_hist_visit;
num_min_hist_visit = 4; % 4 history visits
min_num_visit = num_min_hist_visit + 1; % at least one future visit, this means we will include all eyes which have at least 5 visits
is_parse_glaucoma_data = 1;

%% setup output directory
top_out_dir = sprintf('output_predict_cv2_minvisit%d_avg_rnfl', min_num_visit);
mkdir(top_out_dir);

global fp_log;
str = sprintf('%s/log.txt', top_out_dir);
fp_log = fopen(str, 'wt');

data_parse_output_dir = sprintf('%s/data_stat', top_out_dir);
mkdir(data_parse_output_dir);

%% parse glaucoma Boston dataset

load ../../../onh_registration_clustering/data_new/avg_rnfl_unif/eye_seq_list.mat
load ../../../onh_registration_clustering/data_new/avg_rnfl_unif/bins.mat

for i = 1:length(eye_seq_list)
    
    for v = 1:eye_seq_list{i}.num_visit
        rnfl = eye_seq_list{i}.visit_list{v}.rnfl_ftr;
        eye_seq_list{i}.visit_list{v}.data(2) = rnfl;
        eye_seq_list{i}.visit_data_list(v,2) = rnfl;
    end
    
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
% data_setting.dim_value_range_ls{2} = [110:(-5):30]; % dim2
data_setting.dim_value_range_ls{2} = fliplr(bins);

%% start training CT-HMMs: 
% 10 fold cross validation
num_cv = 2;
cv_idx_list = [1:1:2];

% num_cv = 10;
% cv_idx_list = [1:1:10];

learn_method = 4; % 1: expm, 2: unif, 3: eigen, 4: hybrid (expm in the first iteration, then eigen for the rest)
is_outer_soft = 1; % 1: soft, 2: hard
func_glaucoma_train_cross_validation(top_out_dir, eye_seq_list, learn_method, is_outer_soft, num_cv, cv_idx_list);

%% prediction test using the pre-trained CT-HMM models
pred_out_dir = top_out_dir;
pretrain_CTHMM_dir = top_out_dir;

func_glaucoma_predict_cross_validation(eye_seq_list, pred_out_dir, pretrain_CTHMM_dir, num_cv);
