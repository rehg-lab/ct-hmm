%% load "UCSFFSL"
function [subject_list] = parse_UCSFFreesurfer_datasheet(datasheet_folder)

str = sprintf('%s/UCSFFSL.xlsx', datasheet_folder);
[num,txt,raw] = xlsread(str);
num_row = size(num, 1);
num_col = size(num, 2);

%% field we want to use
field_list = {'RID', 'VISCODE', 'EXAMDATE', 'ST29SV', 'ST88SV', 'ST49CV', 'ST25CV'};
                                          % lefthippo, 
                                          % righthippo,
                                          % ST49CV:leftpostcentral
                                          % ST25CV:leftfrontalpole

for c = 1:num_col
    if(strcmp(txt(1, c), field_list{4}) == 1)
        lefthippo_field = c;
        break;
    end
end

for c = 1:num_col
    if(strcmp(txt(1, c), field_list{5}) == 1)
        righthippo_field = c;
        break;
    end
end

for c = 1:num_col
    if(strcmp(txt(1, c), field_list{6}) == 1)
        leftpostcentral_field = c;
        break;
    end
end

for c = 1:num_col
    if(strcmp(txt(1, c), field_list{7}) == 1)
        leftfrontalpole_field = c;
        break;
    end
end


%% gather data of same subject together
max_ADNI1_line = 3341;
%max_RID = max(num(1:max_ADNI1_line, 1)); % max RID: 
max_RID = 1435;
num_unique_RID = size(unique(num(1:max_ADNI1_line, 1)), 1); % unique RID: 

subject_list = cell(max_RID, 1);
%% init subjectent list
for p = 1:max_RID
    subject_list{p}.num_visit = 0;
end

%% begin fill visiting data
for r = 1:max_ADNI1_line
    
    %% RID
    RID = num(r, 1);
    
    %% fill data    
    k = subject_list{RID}.num_visit + 1;    
    subject_list{RID}.num_visit = k;
    
    subject_list{RID}.visit_list{k}.VISCODE = txt(r+1, 3);
    
    if (strcmp(subject_list{RID}.visit_list{k}.VISCODE{1}, 'sc') == 1)
        subject_list{RID}.visit_list{k}.VISCODE = 'bl';
    end
    
    subject_list{RID}.visit_list{k}.EXAMDATE = txt(r+1, 2);
    subject_list{RID}.visit_list{k}.lefthippo = num(r, lefthippo_field);
    subject_list{RID}.visit_list{k}.righthippo = num(r, righthippo_field);
    subject_list{RID}.visit_list{k}.leftpostcentral = num(r, leftpostcentral_field);
    subject_list{RID}.visit_list{k}.leftfrontalpole = num(r, leftfrontalpole_field);
    
    
    
    
end

subject_list
