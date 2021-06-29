%% decode the glaucoma dataset
% Run the script to get the most probable state path for each eye.
% Specifying a trained model or doing training again in the script, it will give you the decoded path as the figure attached.
% You can modify some input parameters, such as top_out_dir, datasheet_name, min_num_visit, data_setting.type_name_ls,  data_setting.dim_value_range_ls, run_model_training, dwelling_time_draw_range_list, dwelling_time_draw_color_list.
% If you don't want to train the model again, set  run_model_training = 0, and set  top_out_dir to the folder where your model is.

%clear all; 
close all; clc

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

%% add decoding paths
addpath('../../decode/baseline_SSA');
addpath('../../decode');

%% setup output directory
%top_out_dir = '../../../output_glaucoma_decoding_cal2';
%top_out_dir = '../../../output_glaucoma_decoding';
top_out_dir = 'output_decoding_minvisit5';

mkdir(top_out_dir);

%% please put your datasheet in datasheet_folder 
datasheet_folder = 'Glaucoma_Boston_parse';
addpath(datasheet_folder);

%% specifiy the datasheet name here
%datasheet_name = 'yuying2_cal.xlsx'; % this is the datasheet that Igor sent to me in Aug, 2016 (RNFL already calibrated)
datasheet_name = 'Pitt_cthmm_116_6m.xlsx'; 
datasheet_file = sprintf('%s/%s', datasheet_folder, datasheet_name);

global fp_log;
str = sprintf('%s/log.txt', top_out_dir);
fp_log = fopen(str, 'wt');

%% load eye data from the datasheet (here, we include all eyes with at least 5 visit data)
%num_min_hist_visit = 4; % 4 history visits
min_num_visit = 5; % at least one future visit

%% parse glaucoma Boston dataset
data_parse_output_dir = sprintf('%s/data_stat', top_out_dir);
mkdir(data_parse_output_dir);

is_parse_glaucoma_data = 1;

if (is_parse_glaucoma_data == 1)
    %% select patients
    is_extract_g_gs = 1; % extract glaucoma and glaucoma suspect  
    str_Dx = 'g/gs'; % diagnosis code
    dataset_name = 'Boston';
    datasheet_folder = 'Glaucoma_Boston_parse';       
    addpath(datasheet_folder);
    
    %% calibrated data
    %% change data column to be data1, data2
    [all_eye_list, qual_eye_list, num_visit_list, trace_year_list, data1_list, data2_list, baseline_age_list, age_list] = ...
        func_parse_glaucoma_Boston_datasheet_calibrated_general(datasheet_file, data_parse_output_dir, is_extract_g_gs, min_num_visit);
	
    %[all_eye_list, qual_eye_list, num_visit_list, trace_year_list, VFI_list, RNFL_list, baseline_age_list, age_list] = ...
       % func_parse_glaucoma_Boston_datasheet_calibrated(datasheet_file, data_parse_output_dir, is_extract_g_gs, min_num_visit);
    
    %% original data used in NIPS
    %[all_eye_list, qual_eye_list, num_visit_list, trace_year_list, VFI_list, RNFL_list, baseline_age_list, age_list] = ...
     %   func_parse_glaucoma_Boston_datasheet(datasheet_folder, is_extract_g_gs, min_num_visit);    
    
    %% select training group
    age_l = 0; age_u = 120; % age low bound and upper bound, select all ages
    [eye_seq_list] = func_glaucoma_select_train_visit_age(age_l, age_u, qual_eye_list);        
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
state_definition = 2;
if (state_definition == 1)
%% state definition 1
    data_setting.dim_value_range_ls{1} = [100.5 99.5 98 96 93 90 85 80:(-10):20]; % dim1: VFI
    data_setting.dim_value_range_ls{2} = [130:(-5):50 40 30]; % dim2: RNFL
elseif (state_definition == 2)
    %% state definition 2
    data_setting.dim_value_range_ls{1} = [100.5 99.5 98 96 94 92 90 85 80:(-10):20]; % dim1
    data_setting.dim_value_range_ls{2} = [130:(-5):70 60:(-10):30]; % dim2
end

%% cross validation setting
num_cv = 1;
cv_idx_list = [1];

%% start training CT-HMM
global is_draw_learn_Q_mat;
is_draw_learn_Q_mat = 1;
learn_method = 4; % hybrid method (first iter used expm, then uses eigen for the rest)
is_outer_soft = 1; % soft

%% here, you can train or load the pre-trained model for decoding the sequences

run_model_training = 0; %% changed from 1 to 0
if (run_model_training == 1) % train the model
	func_glaucoma_train_cross_validation(top_out_dir, eye_seq_list, learn_method, is_outer_soft, num_cv, cv_idx_list);
else
	%% load the pre-trained CT-HMM model

    % read the number of iteration from num_iter.txt file	
	str = sprintf('%s/CV_1/num_iter.txt', top_out_dir);    
    fp = fopen(str, 'rt');
    num_iter = fscanf(fp, '%d');
    
	%CTHMM_model_dir = 'C:/output_glaucoma_decoding/CV_1/Iter_19'; % original data
	%CTHMM_model_dir = 'C:/exp_output/output_glaucoma_vis_minvisit5_hybrid_statedef2/age_0_110/CV_1/Iter_18';
    CTHMM_model_dir = sprintf('%s/CV_1/Iter_%d', top_out_dir, num_iter); 
	
    CTHMM_learn_load_para(CTHMM_model_dir);
	
end

global state_list;

%% start decoding every sequence
num_obs_seq = size(eye_seq_list, 1);
global is_decode_simulation;
is_decode_simulation = 0;
global is_glaucoma_decoding;
is_glaucoma_decoding = 1;

%% color setting for visualization of state dwelling time
global dwelling_time_draw_range_list;
global dwelling_time_draw_color_list;
dwelling_time_draw_range_list = [0 2 5 100];
dwelling_time_draw_color_list = [255 0 0 ; 255 0 0; 0 255 0; 0 255 0];

%% setup decoding output directory
decode_fig_dir = sprintf('%s/decoded_traj', top_out_dir);
mkdir(decode_fig_dir);

%% start decoding every sequence
run = 1;
if (run == 1)
	for i = 1:num_obs_seq		
		i		
		%% do continuous-time decoding
		eye_seq_list{i}.has_compute_data_emiss_prob = 0;		
		tStart = tic;    
		expect_dur_method = 1; % method for computing expected dwell time for the state: expm
		[decoded_conti_state_seq, decoded_conti_dwelltime_seq, decoded_obs_state_seq] = CTHMM_decode_complete_path_viterbiexpect(eye_seq_list{i}, expect_dur_method);
		tEnd = toc(tStart);

		span_year = (eye_seq_list{i}.visit_time_list(end) - eye_seq_list{i}.visit_time_list(1)) / 12.0;
		str = sprintf('\nSequence decoding time: %d minutes and %f seconds (visit count = %d, time span = %f years)\n\n', floor(tEnd/60),rem(tEnd,60), eye_seq_list{i}.num_visit, span_year);
		CTHMM_print_log(str);    
		
		eye_seq_list{i}.decoded_conti_state_seq = decoded_conti_state_seq;
		eye_seq_list{i}.decoded_conti_dwelltime_seq = decoded_conti_dwelltime_seq;
		eye_seq_list{i}.decoded_obs_state_seq = decoded_obs_state_seq;
		
		%% visualize the continuous decoded  trajectory        
		out_filename = sprintf('%s/decoded_traj_%03d_%d_%s', decode_fig_dir, i, eye_seq_list{i}.ID, eye_seq_list{i}.eye{1});
		CTHMM_vis_2D_decoded_conti_traj(eye_seq_list{i}, decoded_conti_state_seq, decoded_conti_dwelltime_seq, decoded_obs_state_seq, out_filename);
	end
end

%% save the decoded results and state_list
outfile = sprintf('%s/decoded_eye_seq_list', top_out_dir);
save(outfile, 'eye_seq_list');
%load(outfile);

%% output the sequence index that passed each particular state into a log file
str = sprintf('%s/log_decoded_result.txt', top_out_dir);
fp_decode = fopen(str, 'wt');

%% for each state, output the eye index that passed this state
num_state = size(state_list, 1);

for s = 1:num_state

    %% state index and data range
    fprintf(fp_decode, 's: %d, range: [%.1f-%.1f; %.1f-%.1f]\n', s, state_list{s}.range(1,1), state_list{s}.range(1,2), state_list{s}.range(2,1), state_list{s}.range(2,2));
    
    %% the subject index that passes this state
    for i = 1:num_obs_seq                
        num_decoded_state = length(eye_seq_list{i}.decoded_conti_state_seq);
        decoded_state_seq = eye_seq_list{i}.decoded_conti_state_seq;
        for k = 1:num_decoded_state
            if (decoded_state_seq(k) == s) % if the decoded state sequence has state s, print the subject info into the log result                
                fprintf(fp_decode, '%d (ID:%d, eye:%s)\n', i, eye_seq_list{i}.ID, eye_seq_list{i}.eye{1});
                break;
            end
        end        
    end
    fprintf(fp_decode, '\n');
    
end

fclose(fp_log);
fclose(fp_decode);
