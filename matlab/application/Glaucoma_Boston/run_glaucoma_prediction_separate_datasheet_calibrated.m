close all; clc

%% please put your datasheet in datasheet_folder 
datasheet_folder = 'Glaucoma_Boston_parse';
addpath(datasheet_folder);

train_datasheet_name = 'Pitt_cthmm_116_6m.xlsx';
%test_datasheet_name = 'Pitt_cthmm_116_6m.xlsx'; %Aushim_HMM_FGP_Data_VFI.xlsx
test_datasheet_name = 'Aushim_HMM_FGP_Data_VFI.xlsx';
    
%% for training, set the min_num_train_visit_per_eye
num_min_visit_per_eye_for_train = 5;
%% for prediction, set the num of minimun history visit for testing (prediction)
global num_min_hist_visit;
num_min_hist_visit = 1; % number of history visits
num_min_visit_per_eye_for_test = num_min_hist_visit + 1;

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

%% setup output directory
top_out_dir = sprintf('output_pred_sep_datasheet_trainminv%d_testminv%d', num_min_visit_per_eye_for_train, num_min_visit_per_eye_for_test);
mkdir(top_out_dir);
str = sprintf('%s/log.txt', top_out_dir);

global fp_log;
fp_log = fopen(str, 'wt');

global obs_seq_list;

do_CTHMM_training = 1;

if (do_CTHMM_training == 1) % use training datasheet for training
    
    %% specifiy the training datasheet path
    train_datasheet_file = sprintf('%s/%s', datasheet_folder, train_datasheet_name);

    %% load eye data from the datasheet (here, we include all eyes with at least 5 visit data)
    train_data_parse_output_dir = sprintf('%s/train_data_stat', top_out_dir);
    mkdir(train_data_parse_output_dir);

    %% parse training data     
    is_extract_g_gs = 1; % extract glaucoma and glaucoma suspect  
    %str_Dx = 'g/gs'; % diagnosis code
    disp('Start parsing TRAINING glaucoma datasheet...');    

    [all_eye_list, qual_eye_list, num_visit_list, trace_year_list, data1_list, data2_list, baseline_age_list, age_list] = ...
        func_parse_glaucoma_Boston_datasheet_calibrated_general(train_datasheet_file, train_data_parse_output_dir, is_extract_g_gs, num_min_visit_per_eye_for_train);
    
    %% select training group
    age_l = 0; age_u = 120; % age low bound and upper bound of selected ages
    [train_eye_seq_list, train_trace_year_list, train_num_visit_list] = func_glaucoma_select_train_visit_age(age_l, age_u, qual_eye_list);        
    
    save('train_eye_seq_list', 'train_eye_seq_list');
    save('train_trace_year_list', 'train_trace_year_list');
    save('train_num_visit_list', 'train_num_visit_list');
    obs_seq_list = train_eye_seq_list;
    
    num_cv = 1;
    cv_idx_list = [1];   
    learn_method = 4; % 1: expm, 2: unif, 3: eigen, 4: hybrid (expm in the first iteration, then eigen for the rest)
    is_outer_soft = 1; % 1: soft, 2: hard
    %% note: if set only 1 cv, means we train one model for all data
    func_glaucoma_train_cross_validation(top_out_dir, train_eye_seq_list, learn_method, is_outer_soft, num_cv, cv_idx_list);
    pretrain_CTHMM_dir = sprintf('%s/CV_1', top_out_dir);
        
    %% train Bayes prior model here
    num_train_eye = size(obs_seq_list, 1);
    train_idx_list = [1:1:num_train_eye]';
    
    %% now, train bayesian prior model using trainning set
    %% first, compute global linear regression parameters for all seqs (on all visits)
    for s = 1:num_train_eye  
        %% compute global linear regression on all visits
        begin_visit_idx = 1;
        end_visit_idx = obs_seq_list{s}.num_visit;
        time_origin_idx = 1;
        [global_LR_regress, global_LR_err_sigma_vec] = func_compute_linear_regress_para(obs_seq_list{s}, begin_visit_idx, end_visit_idx, time_origin_idx);    
        obs_seq_list{s}.global_LR_regress = global_LR_regress;
        obs_seq_list{s}.global_LR_err_sigma_vec = global_LR_err_sigma_vec;        
    end

    is_bayes_joint_inference = 1;    
    [pretrain_bayes_prior_model, pretrain_bayes_mean_sigma, pretrain_bayes_std_sigma] = func_train_Bayes_prior_model(train_idx_list, is_bayes_joint_inference);
    
else
    
    %% you may also load a pretrained CT-HMM model and a bayes model, by specifying the loading path and variables
    
    pretrain_CTHMM_dir = sprintf('%s/CV_1', top_out_dir);     
    %pretrain_CTHMM_dir = 'output_predict_separate_datasheet_trainminv5_testminv5/CV_1';
    
    %% load a pretrained bayes prior model    
    saved_all_variables = sprintf('%s/saved_all_variables', top_out_dir);
    temp = load(saved_all_variables);    
    pretrain_bayes_prior_model = temp.pretrain_bayes_prior_model;    
    pretrain_bayes_mean_sigma = temp.pretrain_bayes_mean_sigma;
    
end

do_CTHMM_prediction = 1;
if (do_CTHMM_prediction == 1)

    %% prediction test using the pre-trained CT-HMM models
    pred_out_dir = sprintf('%s/pred_result', top_out_dir);
    mkdir(pred_out_dir);
    
    %% parse testing datasheet
    %% specifiy the testing datasheet path
    test_datasheet_file = sprintf('%s/%s', datasheet_folder, test_datasheet_name);

    %% load eye data from the datasheet (here, we include all eyes with at least 5 visit data)
    test_data_parse_output_dir = sprintf('%s/test_data_stat', top_out_dir);
    mkdir(test_data_parse_output_dir);

    %% parse glaucoma Boston dataset
    %% select patients
    is_extract_g_gs = 1; % extract glaucoma and glaucoma suspect      
    disp('Start parsing TESTING glaucoma datasheet...');   
    
    [all_eye_list, qual_eye_list, num_visit_list, trace_year_list, data1_list, data2_list, baseline_age_list, age_list] = ...
        func_parse_glaucoma_Boston_datasheet_calibrated_general(test_datasheet_file, test_data_parse_output_dir, is_extract_g_gs, num_min_visit_per_eye_for_test);
    
    
    %% select testing age group
    age_l = 0; age_u = 120; % age low bound and upper bound, select all ages
    [test_eye_seq_list, test_trace_year_list, test_num_visit_list] = func_glaucoma_select_train_visit_age(age_l, age_u, qual_eye_list);            
    
    %% do prediction for the testing eye list
    func_glaucoma_predict_separate_datasheet(test_eye_seq_list, pred_out_dir, pretrain_CTHMM_dir, pretrain_bayes_prior_model, pretrain_bayes_mean_sigma);
    
end

%% save all variables in the workspace to this variables, you can load later to analyze/use the variables inside
saved_all_variables = sprintf('%s/saved_all_variables', top_out_dir);
save(saved_all_variables);

fclose(fp_log);