close all; clc

%% load eye data from the datasheet (here, we include all eyes with at least 5 visit data)
global num_min_hist_visit;
num_min_hist_visit = 4; % 4 history visits
min_num_visit = num_min_hist_visit + 1; % at least one future visit, this means we will include all eyes which have at least 5 visits

%% setup output directory
top_out_dir = sprintf('output_predict_cv2_minvisit%d_delta_rnfl',min_num_visit);
mkdir(top_out_dir);

global fp_log;
str = sprintf('%s/log.txt', top_out_dir);
fp_log = fopen(str, 'wt');

data_parse_output_dir = sprintf('%s/data_stat', top_out_dir);
mkdir(data_parse_output_dir);

%% compute the delta rnfl dataset

temp = load('./data_orig/Boston_qual_eye_list.mat');
eye_seq_list_orig = temp.qual_eye_list;

eye_seq_list = eye_seq_list_orig;

for i = 1:length(eye_seq_list)
    
    rnfl = eye_seq_list_orig{i}.visit_data_list(:,2);
    rnfl_delta = zeros(size(rnfl));
    
    for v = 1:length(rnfl)
        rnfl_delta(v) = rnfl(1) - rnfl(v);
    end
    
    eye_seq_list{i}.visit_data_list(:,2) = rnfl_delta;
    
    for v = 1:length(eye_seq_list{i}.visit_list)
        eye_seq_list{i}.visit_list{v}.data(2) = eye_seq_list{i}.visit_data_list(v,2);        
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

data_setting.dim_value_range_ls{1} = [100 99 98 96 94 92 90 85 80 75 70 65 60:(-10):0]; % dim1 vfi
% data_setting.dim_value_range_ls{2} = [0:4:40 70 80]; % dim 2 rnfl
data_setting.dim_value_range_ls{2} = [-80:20:0 5:5:50 52];

%% start training CT-HMMs: 
% 10 fold cross validation
% num_cv = 1;
% cv_idx_list = 1:1:1;

num_cv = 2;
cv_idx_list = [1:1:2];

% num_cv = 10;
% cv_idx_list = [1:1:10];


learn_method = 4; % 1: expm, 2: unif, 3: eigen, 4: hybrid (expm in the first iteration, then eigen for the rest)
is_outer_soft = 2; % 1: soft, 2: hard
% func_glaucoma_train_cross_validation(top_out_dir, eye_seq_list, learn_method, is_outer_soft, num_cv, cv_idx_list);

%% prediction test using the pre-trained CT-HMM models
pred_out_dir = top_out_dir;
pretrain_CTHMM_dir = top_out_dir;

func_glaucoma_predict_cross_validation(eye_seq_list, pred_out_dir, pretrain_CTHMM_dir, num_cv);
