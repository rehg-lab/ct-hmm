%function [merged_subject_list normal_subject_list MCIAD_subject_list] = parse_ADNI1_data(datasheet_folder)
function [merged_subject_list] = parse_ADNI1_data(datasheet_folder)

global out_dir;


%% merge patient's data from all four datasheets
%% load "diagnostic summary", "ADASSCORES", "UPENNBIOMK6_V2", "UCSFFSL"

[subject_list1] = parse_diagnosis_datasheet(datasheet_folder);
[subject_list2] = parse_ADASscore_datasheet(datasheet_folder);
[subject_list3] = parse_UPennBiomark_datasheet(datasheet_folder);
[subject_list4] = parse_UCSFFreesurfer_datasheet(datasheet_folder);

max_subject_RID = 1435;
%merged_subject_list = cell(max_subject_RID, 1);
merged_subject_list = subject_list1;

num_subject_complete_data = 0;

for RID = 1:max_subject_RID
    
    %% diagnosis
    num_visit1 = subject_list1{RID}.num_visit;
    for v1 = 1:num_visit1
        merged_subject_list{RID}.visit_list{v1}.has_diagnosis = 1;
        merged_subject_list{RID}.visit_list{v1}.has_cognition = 0;
        merged_subject_list{RID}.visit_list{v1}.has_biochemical = 0;
        merged_subject_list{RID}.visit_list{v1}.has_brainMRI = 0;
        merged_subject_list{RID}.visit_list{v1}.data = [nan nan nan nan nan];
        
        VISCODE = merged_subject_list{RID}.visit_list{v1}.VISCODE{1};
        if (strcmp(VISCODE, 'bl') == 1)
            month = 0;
            merged_subject_list{RID}.visit_list{v1}.month = month;
        elseif (VISCODE(1) == 'm')
            month = sscanf(VISCODE, 'm%d');
            merged_subject_list{RID}.visit_list{v1}.month = month;
        else
            str = sprintf('unknown VISCODE=%s\n', VISCODE);
            merged_subject_list{RID}.visit_list{v1}.month = -1;            
        end                                
    end

    if (num_visit1 == 0)
        continue;
    end
    
    %% ADAS
    num_visit2 = subject_list2{RID}.num_visit;
    for v2 = 1:num_visit2
        for v1 = 1:num_visit1
            %% fill TOTAL11
            if (strcmp(subject_list2{RID}.visit_list{v2}.VISCODE, subject_list1{RID}.visit_list{v1}.VISCODE) == 1)                
                merged_subject_list{RID}.visit_list{v1}.TOTAL11 = subject_list2{RID}.visit_list{v2}.TOTAL11;                
                
                if (merged_subject_list{RID}.visit_list{v1}.TOTAL11 < 0)
                    merged_subject_list{RID}.visit_list{v1}.TOTAL11 = 0;
                end
                
                merged_subject_list{RID}.visit_list{v1}.has_cognition = 1;                
                merged_subject_list{RID}.visit_list{v1}.data(1) = merged_subject_list{RID}.visit_list{v1}.TOTAL11;                
                break;
            end
        end
    end

    %% BIOMARK
    num_visit3 = subject_list3{RID}.num_visit;
    for v3 = 1:num_visit3
        for v1 = 1:num_visit1
            %% fill ABETA and TAU
            if (strcmp(subject_list3{RID}.visit_list{v3}.VISCODE, subject_list1{RID}.visit_list{v1}.VISCODE) == 1)
                merged_subject_list{RID}.visit_list{v1}.TAU = subject_list3{RID}.visit_list{v3}.TAU;
                merged_subject_list{RID}.visit_list{v1}.ABETA = subject_list3{RID}.visit_list{v3}.ABETA;                
                merged_subject_list{RID}.visit_list{v1}.has_biochemical = 1;                                
                merged_subject_list{RID}.visit_list{v1}.data(2) = merged_subject_list{RID}.visit_list{v1}.TAU;
                merged_subject_list{RID}.visit_list{v1}.data(3) = merged_subject_list{RID}.visit_list{v1}.ABETA;                
                break;
            end
        end
    end
    
    %% MRI
    num_visit4 = subject_list4{RID}.num_visit;
    for v4 = 1:num_visit4
        for v1 = 1:num_visit1
            %% fill lefthippo and righthippo
            if (strcmp(subject_list4{RID}.visit_list{v4}.VISCODE, subject_list1{RID}.visit_list{v1}.VISCODE) == 1)
                
                merged_subject_list{RID}.visit_list{v1}.lefthippo = subject_list4{RID}.visit_list{v4}.lefthippo;
                merged_subject_list{RID}.visit_list{v1}.righthippo = subject_list4{RID}.visit_list{v4}.righthippo;
                merged_subject_list{RID}.visit_list{v1}.leftpostcentral = subject_list4{RID}.visit_list{v4}.leftpostcentral;
                merged_subject_list{RID}.visit_list{v1}.leftfrontalpole = subject_list4{RID}.visit_list{v4}.leftfrontalpole;
                
                merged_subject_list{RID}.visit_list{v1}.has_brainMRI = 1;                                                
                merged_subject_list{RID}.visit_list{v1}.data(4) = merged_subject_list{RID}.visit_list{v1}.lefthippo;
                merged_subject_list{RID}.visit_list{v1}.data(5) = merged_subject_list{RID}.visit_list{v1}.righthippo;            
                merged_subject_list{RID}.visit_list{v1}.data(6) = merged_subject_list{RID}.visit_list{v1}.leftpostcentral;            
                merged_subject_list{RID}.visit_list{v1}.data(7) = merged_subject_list{RID}.visit_list{v1}.leftfrontalpole;            
                                
                break;
            end
        end
    end
    
    if (num_visit1 > 0 && num_visit2 > 0 && num_visit3 > 0 && num_visit4 > 0)
        num_subject_complete_data = num_subject_complete_data + 1;
    end

end 

num_subject_complete_data

%% plot histogram for each group and for each data type
data_list = cell(3, 7); % 3 groups, 5 data type

for RID = 1:max_subject_RID    

    if (merged_subject_list{RID}.num_visit > 0)   
        
        for v = 1:merged_subject_list{RID}.num_visit
            
            DXCURREN = merged_subject_list{RID}.visit_list{v}.DXCURREN; % diagnosis, group
            
            if (merged_subject_list{RID}.visit_list{v}.has_cognition == 1)
                TOTAL11 =  merged_subject_list{RID}.visit_list{v}.TOTAL11;
                D = data_list{DXCURREN, 1};
                data_list{DXCURREN, 1} = [D TOTAL11];
            end
            if (merged_subject_list{RID}.visit_list{v}.has_biochemical == 1)
                TAU =  merged_subject_list{RID}.visit_list{v}.TAU;
                D = data_list{DXCURREN, 2};
                data_list{DXCURREN, 2} = [D TAU];
                ABETA =  merged_subject_list{RID}.visit_list{v}.ABETA;
                D = data_list{DXCURREN, 3};
                data_list{DXCURREN, 3} = [D ABETA];
            end
            if (merged_subject_list{RID}.visit_list{v}.has_brainMRI == 1)
                
                lefthippo =  merged_subject_list{RID}.visit_list{v}.lefthippo;
                D = data_list{DXCURREN, 4};
                data_list{DXCURREN, 4} = [D lefthippo];
                
                righthippo =  merged_subject_list{RID}.visit_list{v}.righthippo;
                D = data_list{DXCURREN, 5};
                data_list{DXCURREN, 5} = [D righthippo];
                
                leftpostcentral =  merged_subject_list{RID}.visit_list{v}.leftpostcentral;
                D = data_list{DXCURREN, 6};
                data_list{DXCURREN, 6} = [D leftpostcentral];
                
                leftfrontalpole =  merged_subject_list{RID}.visit_list{v}.leftfrontalpole;
                D = data_list{DXCURREN, 7};
                data_list{DXCURREN, 7} = [D leftfrontalpole];
                
                
            end

        end
    end
end


data_type_ls = {'ADAS-Cog', 'Tau', 'ABeta', 'LeftHippo', 'RightHippo', ...
    'LeftPostcentral', 'LeftFrontalpole'};

for type = 1:7
   
figure,
hist(data_list{1,type});
hold on;
hist(data_list{2,type});
hist(data_list{3,type});
h = findobj(gca,'Type','patch');
display(h)
set(h(1),'FaceColor','r','EdgeColor','k');
set(h(2),'FaceColor','b','EdgeColor','k');
set(h(3),'FaceColor','g','EdgeColor','k');
title(data_type_ls{type});

str = sprintf('%s/%s_hist.png', out_dir, data_type_ls{type});
saveas(gcf, str, 'png');
    
end


%% separate subjects to normal and MCI/AD groups
% normal_subject_list = cell(max_subject_RID, 1);
% num_normal = 0;
% MCIAD_subject_list = cell(max_subject_RID, 1);
% num_MCIAD = 0;
% for RID = 1:max_subject_RID    
%     if (merged_subject_list{RID}.num_visit > 0)
%         if (merged_subject_list{RID}.visit_list{1}.DXCURREN == 1)  % normal
%             num_normal = num_normal + 1;
%             normal_subject_list{num_normal} = merged_subject_list{RID};
%             normal_subject_list{num_normal}.RID = RID;
%         else % MCI or AD
%             num_MCIAD = num_MCIAD + 1;
%             MCIAD_subject_list{num_MCIAD} = merged_subject_list{RID};
%             MCIAD_subject_list{num_MCIAD}.RID = RID;
%         end
%     end
% end


%% check that the visit month is increasing
% for s = 1:num_include_subject    
%     num_visit = include_subject_list{s}.num_visit;
%     pre_month = include_subject_list{s}.visit_list{1}.month;     
%     for v = 2:num_visit
%         cur_month = include_subject_list{s}.visit_list{v}.month;
%         if (cur_month < pre_month)  % need swapping
%            disp('month ordering problem');
%         end
%         pre_month = cur_month; 
%     end
% end