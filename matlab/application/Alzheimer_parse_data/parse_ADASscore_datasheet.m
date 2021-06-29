%% load "ADASSCORES"
function [subject_list] = parse_ADASscore_datasheet(datasheet_folder)

str = sprintf('%s/ADASSCORES.xlsx', datasheet_folder);
[num,txt,raw] = xlsread(str);
num_row = size(num, 1);
num_col = size(num, 2);

%% field we want to use
field_list = {'RID', 'VISCODE', 'EXAMDATE', 'TOTAL11'};

for c = 1:num_col
    if(strcmp(txt(1, c), field_list{4}) == 1)
        total11_field = c;
        break;
    end
end


%% gather data of same subject together
max_ADNI1_line = num_row;
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
    
    subject_list{RID}.visit_list{k}.VISCODE = txt(r+1, 4);
    subject_list{RID}.visit_list{k}.EXAMDATE = txt(r+1, 7);
    subject_list{RID}.visit_list{k}.TOTAL11 = num(r, total11_field);

end

subject_list
