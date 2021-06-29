function func_write_excel_eye_list(qual_eye_list, out_xlsx_filename)


%% write the information of each included eye in an excel file
column_names = {'subject', 'eye', 'baseline_age', 'num of visits', 'yrsfrombaseline_list', 'Dx_list', 'data1_list', 'data2_list'};
num_column = 8;
num_qual_eye = size(qual_eye_list, 1);
A = cell(num_qual_eye + 1, num_column);

A(1, :) = column_names;

for  s = 1:num_qual_eye

    row = s+1;
    
    % subject ID
    A{row, 1} = qual_eye_list{s}.ID;
    % eye
    A{row, 2} = qual_eye_list{s}.eye{1};
    % baseline age
    A{row, 3} = qual_eye_list{s}.visit_list{1}.age;
    % num of included visits
    A{row, 4} = qual_eye_list{s}.num_visit;
    
    str_yrsfrombaseline_list = '';
    str_Dx_list = '';
    str_data1_list = '';
    str_data2_list = '';
    
    for v = 1:qual_eye_list{s}.num_visit;
       % yrsfrombaseline_list
       temp = sprintf('%.1f,', (qual_eye_list{s}.visit_list{v}.day - qual_eye_list{s}.visit_list{1}.day) / 365.0);
       str_yrsfrombaseline_list = strcat(str_yrsfrombaseline_list, temp);
       % Dx_list
       temp = sprintf('%s,', qual_eye_list{s}.visit_list{v}.Dx{1});
       str_Dx_list = strcat(str_Dx_list, temp);
       % data1 list
       temp = sprintf('%f,', qual_eye_list{s}.visit_list{v}.data1);
       str_data1_list = strcat(str_data1_list, temp);
       % data2 list       
       temp = sprintf('%f,', qual_eye_list{s}.visit_list{v}.data2);
       str_data2_list = strcat(str_data2_list, temp);
    end
    A{row, 5} = str_yrsfrombaseline_list;
    A{row, 6} = str_Dx_list;
    A{row, 7} = str_data1_list;
    A{row, 8} = str_data2_list;        
end

%out_xlsx_filename = sprintf('%s/selected_subjects.xlsx', output_dir);
% xlswrite(out_xlsx_filename,A); %% changed by supriya
save('out_xlsx_filename','A')
