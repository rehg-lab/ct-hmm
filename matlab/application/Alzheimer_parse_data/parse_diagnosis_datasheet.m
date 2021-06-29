%% load "diagnostic summary"
function [subject_list] = parse_diagnosis_datasheet(datasheet_folder)


str = sprintf('%s/DXSUM_PDXCONV_ADNIALL.xlsx', datasheet_folder);
[num,txt,raw] = xlsread(str);
num_row = size(num, 1);
num_col = size(num, 2);

%% field we want to use
field_list = {'Phase', 'RID', 'VISCODE', 'EXAMDATE', 'DXCURREN'};

%% gather data of same subject together
max_ADNI1_line = 3856;
max_RID = max(num(1:max_ADNI1_line, 2)); % max RID: 1435
num_unique_RID = size(unique(num(1:max_ADNI1_line, 2)), 1); % unique RID : 819

subject_list = cell(max_RID, 1);
%% init subjectent list
for p = 1:max_RID
    subject_list{p}.num_visit = 0;    
end

%% begin fill visiting data
for r = 1:max_ADNI1_line
    
    %% RID
    RID = num(r, 2);
    
    %% fill data    
    k = subject_list{RID}.num_visit + 1;    
    subject_list{RID}.num_visit = k;
    subject_list{RID}.visit_list{k}.VISCODE = txt(r+1, 5);
    subject_list{RID}.visit_list{k}.EXAMDATE = txt(r+1, 9);
    subject_list{RID}.visit_list{k}.DXCURREN = num(r, 10);

end

subject_list
