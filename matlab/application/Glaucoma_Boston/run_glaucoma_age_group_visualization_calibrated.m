close all; clc

%% please put your datasheet in datasheet_folder 
datasheet_folder = 'Glaucoma_Boston_parse';
addpath(datasheet_folder);

%% specifiy the datasheet name here
%datasheet_name = 'yuying2_cal.xlsx'; % this is the datasheet that Igor sent to me in Aug, 2016 (RNFL already calibrated)
datasheet_name = 'Pitt_cthmm_116_6m.xlsx';
datasheet_file = sprintf('%s/%s', datasheet_folder, datasheet_name);

%% load eye data from the datasheet
min_num_visit = 5;   % minimum number of visits for each eye to be included in the model training and visualization
is_parse_glaucoma_data = 1;

%% setup your output directory
%top_out_dir = sprintf('output_vis_minvisit%d_PittNew_tol10-6_def7_0423', min_num_visit);
%top_out_dir = sprintf('output_vis_minvisit%d_PittNew_tol10-6_def4_0427_2018_hybrid', min_num_visit);
top_out_dir = sprintf('output_vis_minvisit%d_PittNew_tol10-6_def4_0427_2018_expm_rerun', min_num_visit);
mkdir(top_out_dir);

global fp_log;
log_filename = sprintf('%s/log.txt', top_out_dir);
fp_log = fopen(log_filename, 'wt');

%% parse glaucoma Boston datasheet
data_parse_output_dir = sprintf('%s/data_stat', top_out_dir); % output data statistics
mkdir(data_parse_output_dir);

if (is_parse_glaucoma_data == 1)
	is_extract_g_gs = 1;   % extract glaucoma and glaucoma suspect
    str_Dx = 'g/gs'; % diagnosis code
    
    disp('Start parsing glaucoma datasheet...');
    
	%% here, parse the datasheet and output the data statistics
    %[all_eye_list, qual_eye_list, num_visit_list, trace_year_list, VFI_list, RNFL_list, baseline_age_list, age_list] = ...
    %    func_parse_glaucoma_Boston_datasheet_calibrated(datasheet_file, data_parse_output_dir, is_extract_g_gs, min_num_visit);
	
	[all_eye_list, qual_eye_list, num_visit_list, trace_year_list, data1_list, data2_list, baseline_age_list, age_list] = ...
        func_parse_glaucoma_Boston_datasheet_calibrated_general(datasheet_file, data_parse_output_dir, is_extract_g_gs, min_num_visit);
		
				
    save('qual_eye_list', 'qual_eye_list');	
else % load the parsed data
    temp = load('qual_eye_list');
    qual_eye_list = temp.qual_eye_list;
end 

%% setup the CT-HMM learning algorithm
% learn_method = 4; % % 1: expm, 2: unif, 3: eigen, 4: hybrid (expm in the first iteration, then eigen for the rest)

learn_method = 1; % % 1: expm, 2: unif, 3: eigen, 4: hybrid (expm in the first iteration, then eigen for the rest)

is_outer_soft = 1; % 1: soft, 2: hard method


%% data setting
global data_setting;
data_setting.dim = 2;
data_setting.type_name_ls = {'VFI', 'RNFL'};

%% setup state grid
state_definition = 4;

%% you can setup your state definition (grid) here
if (state_definition == 1)
	%% state definition 1
    data_setting.dim_value_range_ls{1} = [100.5 99.5 98 96 93 90 85 80:(-10):20]; % dim1: VFI
    data_setting.dim_value_range_ls{2} = [130:(-5):50 40 30]; % dim2: RNFL
elseif (state_definition == 2)
    %% state definition 2
    %data_setting.dim_value_range_ls{1} = [100.5 99.5 98 96 94 92 90 85 80:(-10):20]; % dim1
    %data_setting.dim_value_range_ls{2} = [130:(-5):70 60:(-10):30]; % dim2
    data_setting.dim_value_range_ls{1} = [100.5 99.5 98 96 94 92 90 85 80:(-10):0]; % dim1
    data_setting.dim_value_range_ls{2} = [130:(-5):70 60:(-10):30]; % dim2
elseif (state_definition == 3) % added on 4.27.2018
	data_setting.dim_value_range_ls{1} = [100.5 99.5 98.5 97.5 96.5 94.5 93.5 92.5 91.5 90.5 88:(-2):70 60:(-10):0]; % dim1    
    data_setting.dim_value_range_ls{2} = [110:(-5):30]; % dim2 
elseif (state_definition == 4) % added on 4.27.2018
	data_setting.dim_value_range_ls{1} = [100.5 99.5 98.5 97.5 96.5 94.5 93.5 92.5 91.5 90.5 88:(-2):70 60:(-10):0]; % dim1    
    data_setting.dim_value_range_ls{2} = [110:(-5):100 97.5:(-2.5):50 45:(-5):30]; % dim2    
end

%definition 7
%data_setting.dim_value_range_ls{1} = [100.5 99.5 98 96 94 92 90 85 80 75 70:(-10):0]; % dim1
%data_setting.dim_value_range_ls{2} = [110:(-5):30]; % dim2    

%% for visualization, draw data origin at the left top corner
data_setting.draw_origin_lefttop = 1;

%% this means just train 1 model for all the data (no cross validation)
num_cv = 1; % this means just train 1 model
cv_idx_list = 1;

%% select training group
%% all age

age_group_list = [0 120; 0 70; 70 120]; % age group [0 120], [0 70], [70 120]

for age_gr_idx = 1:1:1 % run for just first age group
%for age_group = 1:1:3  % run for all three specified age groups

	age_l = age_group_list(age_gr_idx, 1);
	age_u = age_group_list(age_gr_idx, 2);
	
	%% select visit data that are within the specified age range
	[select_eye_seq_list] = func_glaucoma_select_train_visit_age(age_l, age_u, qual_eye_list);
	age_top_out_dir = sprintf('%s/age_%d_%d', top_out_dir, age_l, age_u);
	mkdir(age_top_out_dir);
	
	out_xlsx_filename = sprintf('%s/selected_subjects.xlsx', age_top_out_dir);
	func_write_excel_eye_list(select_eye_seq_list, out_xlsx_filename);
	func_glaucoma_train_cross_validation(age_top_out_dir, select_eye_seq_list, learn_method, is_outer_soft, num_cv, cv_idx_list);

end

fclose(fp_log);
