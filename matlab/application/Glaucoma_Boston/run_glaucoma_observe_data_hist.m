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

%% setup your output directory
top_out_dir = sprintf('data_observe');
mkdir(top_out_dir);

global fp_log;

log_filename = sprintf('%s/log.txt', top_out_dir);
fp_log = fopen(log_filename, 'wt');

%% parse glaucoma Boston datasheet
data_parse_output_dir = top_out_dir;

is_extract_g_gs = 1;   % extract glaucoma and glaucoma suspect
str_Dx = 'g/gs'; % diagnosis code

disp('Start parsing glaucoma datasheet...');
    
%% here, parse the datasheet and output the data statistics
[all_eye_list, qual_eye_list, num_visit_list, trace_year_list, data1_list, data2_list, baseline_age_list, age_list] = ...
        func_parse_glaucoma_Boston_datasheet_calibrated_general(datasheet_file, data_parse_output_dir, is_extract_g_gs, min_num_visit);
		

data_setting.type_name_ls = {'VFI', 'RNFL'};

%% setup state grid
data_setting.dim_value_range_ls{1} = [100.5 99.5 98 96 94 92 90 85 80 75 70:(-10):0]; % dim1: VFI
data_setting.dim_value_range_ls{2} = [110:(-5):30]; % dim2: RNFL   

dim1_edges = data_setting.dim_value_range_ls{1};
dim2_edges = data_setting.dim_value_range_ls{2};
num_dim1_bin = length(dim1_edges)-1;
num_dim2_bin = length(dim2_edges)-1;

raw_data_2D_hist = zeros(num_dim1_bin, num_dim2_bin);
num_raw_data = length(data1_list);

for i = 1:num_raw_data
    % find which bin the data is in
    for x = 1:num_dim1_bin        
        if ((data1_list(i) <= dim1_edges(x)) && (data1_list(i) > dim1_edges(x+1)))         
            for y = 1:num_dim2_bin                
                if ((data2_list(i) <= dim2_edges(y)) && (data2_list(i) > dim2_edges(y+1)))         
                    % find the bin
                    raw_data_2D_hist(x, y) = raw_data_2D_hist(x, y) + 1; 
                    break;                    
                end
            end        
        end        
    end    
end

func_vis_2D_raw_data_hist(top_out_dir, dim1_edges, dim2_edges, raw_data_2D_hist);

%histogram2(data1_list,data2_list,data_setting.dim_value_range_ls{1},data_setting.dim_value_range_ls{2}, 'DisplayStyle','tile','ShowEmptyBins','on');
%colorbar

