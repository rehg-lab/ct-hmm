function [eye_list, subject_list, num_visit_list, trace_year_list, VFI_list, RNFL_list, baseline_age_list, age_list] = ...    
parse_glaucoma_Boston_datasheet(datasheet_folder)

str = sprintf('%s\\longitudinal dataset new.xlsx',datasheet_folder);

[num,txt,raw] = xlsread(str);
num_row = size(num, 1);
num_col = size(num, 2);

%% field we want to use
field_list = {'subject', 'eye', 'daysfrombaseline', 'Dx', ...
              'vf', 'VFI', ...
              'oct_prototype', 'oct_1_2', 'oct_stratus', 'oct_cirrus', ...
              'mean_rnfl_prototype', 'mean_rnfl_oct_1_2', 'mean_rnfl_stratus', 'mean_rnfl_cirrus', ...
              'ss_stratus', 'ss_cirrus', ...
              'baseline_age', 'age'};
              
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
max_line = num_row;
max_ID = max(num(1:max_line, 1));
num_unique_ID = size(unique(num(1:max_line, 1)), 1);

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

code_count_ls = zeros(16, 1);

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
    
    eye_list{num_eye}.visit_list{k}.month = double(eye_list{num_eye}.visit_list{k}.day) / 30.0;
    
    oct_prototype = num(r, field_colidx_list(7));
    oct_1_2 = num(r, field_colidx_list(8));
    oct_stratus = num(r, field_colidx_list(9));
    oct_cirrus = num(r, field_colidx_list(10));
    mean_rnfl_prototype = num(r, field_colidx_list(11));
    mean_rnfl_oct_1_2 = num(r, field_colidx_list(12));
    mean_rnfl_stratus = num(r, field_colidx_list(13));
    mean_rnfl_cirrus = num(r, field_colidx_list(14));
    ss_stratus = num(r, field_colidx_list(15));
    ss_cirrus = num(r, field_colidx_list(16));
    baseline_age = num(r, field_colidx_list(17));
    age = num(r, field_colidx_list(18));
    
    eye_list{num_eye}.visit_list{k}.oct_prototype = oct_prototype;
    eye_list{num_eye}.visit_list{k}.oct_1_2 = oct_1_2;
    eye_list{num_eye}.visit_list{k}.oct_stratus = oct_stratus;
    eye_list{num_eye}.visit_list{k}.oct_cirrus = oct_cirrus;    
    eye_list{num_eye}.visit_list{k}.mean_rnfl_prototype = mean_rnfl_prototype;
    eye_list{num_eye}.visit_list{k}.mean_rnfl_oct_1_2 = mean_rnfl_oct_1_2;
    eye_list{num_eye}.visit_list{k}.mean_rnfl_stratus = mean_rnfl_stratus;
    eye_list{num_eye}.visit_list{k}.mean_rnfl_cirrus = mean_rnfl_cirrus;
    eye_list{num_eye}.visit_list{k}.ss_stratus = ss_stratus;
    eye_list{num_eye}.visit_list{k}.ss_cirrus = ss_cirrus;
    
    eye_list{num_eye}.visit_list{k}.baseline_age = baseline_age;
    eye_list{num_eye}.visit_list{k}.age = age;
        
    if (oct_prototype ~= 1)
        oct_prototype = 0;
    end
    if (oct_1_2 ~= 1)
        oct_1_2 = 0;
    end
    if (oct_stratus ~= 1)
        oct_stratus = 0;
    end
    if (oct_cirrus ~= 1)
        oct_cirrus = 0;
    end
    
    %% accumulate code
    code = (oct_prototype * 8 + oct_1_2 * 4 + oct_stratus * 2 + oct_cirrus) + 1;
    code_count_ls(code) = code_count_ls(code) + 1;
 
    %% update
    pre_subject = subject;
    pre_eye = eye;
    
end % r

% %% no oct
% str = sprintf('No qual OCT measures: #%d\n', code_count_ls(0+1))
% %% single machine
% str = sprintf('oct_prototype only: #%d\n', code_count_ls(8+1))
% str = sprintf('oct_1_2 only: #%d\n', code_count_ls(4+1))
% str = sprintf('oct_stratus only: #%d\n', code_count_ls(2+1))
% str = sprintf('oct_cirrus only: #%d\n', code_count_ls(1+1))
% %% two machines
% str = sprintf('oct_prototype + oct_1_2: #%d\n', code_count_ls(12+1))
% str = sprintf('oct_prototype + oct_stratus only: #%d\n', code_count_ls(10+1))
% str = sprintf('oct_prototype + oct_cirrus only: #%d\n', code_count_ls(9+1))
% str = sprintf('oct_1_2 + oct_stratus only: #%d\n', code_count_ls(6+1))
% str = sprintf('oct_1_2 + oct_cirrus only: #%d\n', code_count_ls(5+1))
% str = sprintf('oct_stratus + oct_cirrus only: #%d\n', code_count_ls(3+1))
% code_count_ls

%% check all visits from all eyes
for e = 1:num_eye
    num_visit = eye_list{e}.num_visit;
    
    for k = 1:num_visit
        visit = eye_list{e}.visit_list{k};
        
        if (visit.oct_prototype == 1 || visit.oct_1_2 == 1 || visit.oct_stratus == 1 || visit.oct_cirrus == 1)
            %% compute normalized stratus measure
            eye_list{e}.visit_list{k}.has_qual_OCT = 1;
        else
            eye_list{e}.visit_list{k}.has_qual_OCT = 0;            
        end
        
        %% transform to stratus measure
        if (eye_list{e}.visit_list{k}.has_qual_OCT == 1)            
            %S = -12.689 + 1.018 P
            %S = -4.119 + 0.900 O
            %S = -15.149 + 1.242 C
            %% =================================
            if (visit.oct_stratus == 1)
                S = visit.mean_rnfl_stratus;
            elseif (visit.oct_cirrus == 1)
                C = visit.mean_rnfl_cirrus;
                if (isnan(visit.ss_cirrus) == 1) % no ss, use simple equation
                    S = -15.149 + 1.242 * C;
                else % use enhanced equation
                    ss_cirrus = visit.ss_cirrus;
                    ss_stratus = 10;
                    S = -9.761 + (1.778 * ss_stratus - 3.850 * ss_cirrus) + 1.391 * C; 
                end
            elseif (visit.oct_stratus == 1)
                S = visit.mean_rnfl_stratus;
            elseif (visit.oct_1_2 == 1)
                O = visit.mean_rnfl_oct_1_2;
                S = -4.119 + 0.900 * O;
            elseif (visit.oct_prototype == 1)
                P = visit.mean_rnfl_prototype;
                S = -12.689 + 1.018 * P;
            end
            visit.norm_rnfl_stratus = S;
            %% ==================================            
        end
            
        if (eye_list{e}.visit_list{k}.has_qual_OCT == 1 && visit.vf == 1)                
            data1 = visit.VFI;
            data2 = visit.norm_rnfl_stratus;
            eye_list{e}.visit_list{k}.data = [data1 data2];
            eye_list{e}.visit_list{k}.include = 1;
        else
            eye_list{e}.visit_list{k}.include = 0;
        end
        
    end    
end

%% extract qualified data into subject_list
[subject_list, num_visit_list, trace_year_list, VFI_list, RNFL_list, baseline_age_list, age_list] = ...
utility_extract_qual_eye_data(eye_list, fp_log);

%% plot dataset statistics
output_prefix = datasheet_name;
utility_plot_datasheet_statistics(fp_log, output_dir, output_prefix, num_visit_list, trace_year_list, VFI_list, RNFL_list, baseline_age_list, age_list);

fclose(fp_log);

