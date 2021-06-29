function [subject_list] = parse_UPennBiomark_datasheet(datasheet_folder)

%% load "UPENNBIOMK6_V2"
str = sprintf('%s/UPENNBIOMK6_V2.xlsx', datasheet_folder);
[num,txt,raw] = xlsread(str);


%% field we want to use
field_list = {'RID', 'VISCODE', 'EXAMDATE', 'ABETA', 'TAU'};

%% gather data of same subject together
max_ADNI1_line = 324;
max_RID = 1435;

subject_list = cell(max_RID, 1);
%% init subjectent list
for p = 1:max_RID    
    subject_list{p}.num_visit = 0;
    subject_list{p}.data_from_v6 = 0;
    subject_list{p}.data_from_v4 = 0;
    subject_list{p}.data_from_v2 = 0;
    subject_list{p}.data_from_v1 = 0;
end

%% begin fill visiting data
BIOMARK6_RID = [];
for r = 1:max_ADNI1_line
    
    %% RID
    RID = num(r, 1);
    BIOMARK6_RID = [BIOMARK6_RID RID];
    
    %% fill data    
    if (isnan(num(r, 5))==0 && isnan(num(r, 6))==0)
        k = subject_list{RID}.num_visit + 1;    
        subject_list{RID}.num_visit = k;
        subject_list{RID}.data_from_v6 = 1;

        subject_list{RID}.visit_list{k}.VISCODE = txt(r+1, 2);
        subject_list{RID}.visit_list{k}.EXAMDATE = txt(r+1, 3);
        subject_list{RID}.visit_list{k}.ABETA = num(r, 5);
        subject_list{RID}.visit_list{k}.TAU = num(r, 6);
    end

end
num_BIOMARK6_RID = size(unique(BIOMARK6_RID), 2)

%%----------------------------------------------------------
%% load "UPENNBIOMK4"
str = sprintf('%s/UPENNBIOMK4.csv', datasheet_folder);
[num,txt,raw] = xlsread(str);
max_ADNI1_line = size(num, 1);
BIOMARK4_RID = [];
for r = 1:max_ADNI1_line
    
    %% RID
    RID = num(r, 1);       
     
    if (subject_list{RID}.data_from_v6 == 0 && isnan(num(r, 5))==0 && isnan(num(r, 4))==0)
        
        BIOMARK4_RID = [BIOMARK4_RID RID];
        k = subject_list{RID}.num_visit + 1;
        subject_list{RID}.num_visit = k;
        subject_list{RID}.data_from_v4 = 1;
        subject_list{RID}.visit_list{k}.VISCODE = txt(r+1, 2);
        subject_list{RID}.visit_list{k}.EXAMDATE = txt(r+1, 3);
        
        subject_list{RID}.visit_list{k}.ABETA = num(r, 5);
        subject_list{RID}.visit_list{k}.TAU = num(r, 4);
        
    end
end
num_BIOMARK4_RID = size(unique(BIOMARK4_RID), 2)

%%----------------------------------------------------------
%% load "UPENNBIOMK2"
str = sprintf('%s/UPENNBIOMK2.csv', datasheet_folder);
[num,txt,raw] = xlsread(str);
max_ADNI1_line = size(num, 1);
BIOMARK2_RID = [];

for r = 1:max_ADNI1_line
    
    %% RID
    RID = num(r, 1);   
    
    if (subject_list{RID}.data_from_v6 == 0 && subject_list{RID}.data_from_v4 == 0 && isnan(num(r, 5))==0 && isnan(num(r, 6))==0)
        
        BIOMARK2_RID = [BIOMARK2_RID RID];
        
        k = subject_list{RID}.num_visit + 1;
        subject_list{RID}.num_visit = k;
        subject_list{RID}.data_from_v2 = 1;
        subject_list{RID}.visit_list{k}.VISCODE = txt(r+1, 3);
        %subject_list{RID}.visit_list{k}.EXAMDATE = txt(r+1, 3);
        subject_list{RID}.visit_list{k}.ABETA = num(r, 6);
        subject_list{RID}.visit_list{k}.TAU = num(r, 5);
    end
end

num_BIOMARK2_RID = size(unique(BIOMARK2_RID), 2)

%%----------------------------------------------------------
%% load "UPENNBIOMK2"
str = sprintf('%s/UPENNBIOMK.csv', datasheet_folder);
[num,txt,raw] = xlsread(str);
max_ADNI1_line = size(num, 1);
BIOMARK1_RID = [];

for r = 1:max_ADNI1_line
    
    %% RID
    RID = num(r, 2);   
    if (subject_list{RID}.data_from_v6 == 0 && subject_list{RID}.data_from_v4 == 0 && subject_list{RID}.data_from_v2 == 0 && isnan(num(r, 7))==0 && isnan(num(r, 8))==0)
        
        BIOMARK1_RID = [BIOMARK1_RID RID];
        
        k = subject_list{RID}.num_visit + 1;
        subject_list{RID}.num_visit = k;
        subject_list{RID}.data_from_v1 = 1;
        subject_list{RID}.visit_list{k}.VISCODE = txt(r+1, 4);
        %subject_list{RID}.visit_list{k}.EXAMDATE = txt(r+1, 3);
        subject_list{RID}.visit_list{k}.ABETA = num(r, 8);
        subject_list{RID}.visit_list{k}.TAU = num(r, 7);
    end
end

num_BIOMARK1_RID = size(unique(BIOMARK1_RID), 2)

