function [eye_list, qual_eye_list, num_visit_list, trace_year_list, VFI_list, RNFL_list, baseline_age_list, age_list] = ...    
                                                func_parse_glaucoma_Boston_datasheet_calibrated(datasheet_file, output_dir, is_extract_g_gs, min_num_visit)

global fp_log;
                                            
datasheet_name = 'Boston';

if (is_extract_g_gs == 1)
    str_Dx = 'ggs';
else
    str_Dx = 'n';
end

%str = sprintf('%s/longitudinal dataset new.xlsx',datasheet_folder);
%str = sprintf('%s/yuying2_cal.xlsx',datasheet_folder);


[num,txt,raw] = xlsread(datasheet_file);

num_row = size(num, 1);
num_col = size(num, 2);

%% field we want to use
field_list = {'subject', 'eye', 'daysfrombaseline', 'Dx', 'vf', 'VFI', 'rnfl', 'baseline_age', 'age'};              
              %'oct_prototype', 'oct_1_2', 'oct_stratus', 'oct_cirrus', 
              %'mean_rnfl_prototype', 'mean_rnfl_oct_1_2', 'mean_rnfl_stratus', 'mean_rnfl_cirrus',
              %'ss_stratus', 'ss_cirrus',               
              
num_used_field = size(field_list, 2);
field_colidx_list = zeros(num_used_field, 1);

%% find field position
for f = 1:num_used_field
    for c = 1:num_col
        if (strcmp(txt(1, c), field_list{f})==1)
            field_colidx_list(f) = c;
            break;
        end
    end
end
          
%% gather data of same eye together
subject_field_idx = field_colidx_list(1);

max_line = num_row;
max_ID = max(num(1:max_line, subject_field_idx));
num_unique_ID = size(unique(num(1:max_line, subject_field_idx)), 1);

%% init subject list
eye_list = cell(num_unique_ID*2, 1);

%% init eye list
for e = 1:num_unique_ID*2
    eye_list{e}.num_visit = 0;
end

%% begin fill visiting data
num_eye = 0;

pre_subject = -1;
pre_eye = 'x';

for r = 1:max_line    
    
    %% subject ID
    subject = num(r, field_colidx_list(1));
    
    %% eye
    eye = txt(r+1, field_colidx_list(2));
    
    if (pre_subject ~= subject || strcmp(pre_eye, eye) == 0)  % a new eye
        num_eye = num_eye + 1;
        eye_list{num_eye}.ID = subject;
        eye_list{num_eye}.eye = eye;
    end

    %% fill data
    k = eye_list{num_eye}.num_visit + 1;
    eye_list{num_eye}.num_visit = k;
                           
    eye_list{num_eye}.visit_list{k}.day = num(r, field_colidx_list(3));           
    eye_list{num_eye}.visit_list{k}.Dx = txt(r+1, field_colidx_list(4));  % we want g or gs
    eye_list{num_eye}.visit_list{k}.vf = num(r, field_colidx_list(5));
    eye_list{num_eye}.visit_list{k}.VFI = num(r, field_colidx_list(6));
    
    eye_list{num_eye}.visit_list{k}.month = round(double(eye_list{num_eye}.visit_list{k}.day) / 30.0);    
   
    %eye_list{num_eye}.visit_list{k}.time = double(eye_list{num_eye}.visit_list{k}.day) / 365.0;
    eye_list{num_eye}.visit_list{k}.time = eye_list{num_eye}.visit_list{k}.month;
    
    RNFL = num(r, field_colidx_list(7)); % new 
    
    baseline_age = num(r, field_colidx_list(8));
    age = num(r, field_colidx_list(9));
    
    eye_list{num_eye}.visit_list{k}.RNFL = RNFL; % new
    
    eye_list{num_eye}.visit_list{k}.baseline_age = baseline_age;
    eye_list{num_eye}.visit_list{k}.age = age;
 
    %% update
    pre_subject = subject;
    pre_eye = eye;
    
end % r

%% check all visits from all eyes
for e = 1:num_eye
    
    num_visit = eye_list{e}.num_visit;
    
    for k = 1:num_visit
        
        visit = eye_list{e}.visit_list{k};
                   
        if (visit.vf == 1 && isnan(visit.VFI) == 0 && isnan(visit.RNFL) == 0)                
            data1 = visit.VFI;
            data2 = visit.RNFL;
            eye_list{e}.visit_list{k}.data = [data1 data2];
            eye_list{e}.visit_list{k}.include = 1;
        else
            eye_list{e}.visit_list{k}.include = 0;
        end
        
    end    
end

%str = sprintf('%s/log_statistics.txt', output_dir);
%fp_log = fopen(str, 'wt');

%% extract qualified data into subject_list
[qual_eye_list, num_visit_list, trace_year_list, VFI_list, RNFL_list, baseline_age_list, age_list] = ...
func_extract_qual_eye_data(eye_list, is_extract_g_gs, min_num_visit, fp_log);

%% write the information of each included eye in an excel file
out_xlsx_filename = sprintf('%s/selected_subjects_%s_minvisit%d.xlsx', output_dir, str_Dx, min_num_visit);
func_write_excel_eye_list(qual_eye_list, out_xlsx_filename);

%% plot dataset statistics
output_prefix = sprintf('%s_%s', datasheet_name, str_Dx);
func_plot_datasheet_statistics(output_dir, output_prefix, num_visit_list, trace_year_list, VFI_list, RNFL_list, baseline_age_list, age_list);

