function [include_subject_list] = func_alzheimer_select_ND_subject_data(complete_subject_list, use_cognition, use_brainMRI, useBiochemical, data_idx)

max_subject_RID = size(complete_subject_list, 1);
include_subject_list = cell(max_subject_RID, 1);
num_include_subject = 0;
total_num_include_visit = 0;

dim = size(data_idx, 2);

visit_num_list = [];
trace_year_list = [];

for RID = 1:max_subject_RID
     
     num_visit = complete_subject_list{RID}.num_visit;
     num_include_visit = 0;     
     include_visit_list = [];
     
     for v = 1:num_visit                  
         if (use_cognition == 1 && complete_subject_list{RID}.visit_list{v}.has_cognition == 0)
            continue;
         end
         if (use_brainMRI == 1 && complete_subject_list{RID}.visit_list{v}.has_brainMRI == 0)
            continue;
         end
         if (useBiochemical == 1 && complete_subject_list{RID}.visit_list{v}.has_biochemical == 0)
            continue;
         end
         if (complete_subject_list{RID}.visit_list{v}.DXCURREN == 1) % normal
             continue;
         end
         if (complete_subject_list{RID}.visit_list{v}.month  < 0) % invalid
             continue;
         end
         num_include_visit = num_include_visit + 1;                  
         include_visit_list = [include_visit_list v];
     end
         
     %% add this subject
     %% select at least 1 visits
     if (num_include_visit >= 1)
         
         num_include_subject = num_include_subject + 1;
         
         include_subject_list{num_include_subject}.num_visit = num_include_visit;
         include_subject_list{num_include_subject}.visit_list = cell(num_include_visit,1);
         include_subject_list{num_include_subject}.RID = RID;
         
         include_subject_list{num_include_subject}.visit_time_list = zeros(num_include_visit, 1);
         include_subject_list{num_include_subject}.visit_data_list = zeros(num_include_visit, dim);
                  
         for t = 1:num_include_visit
             v = include_visit_list(t);
             include_subject_list{num_include_subject}.visit_list{t} = complete_subject_list{RID}.visit_list{v};             
         end
         
         trace_year = include_subject_list{num_include_subject}.visit_list{end}.month / 12.0;
         
         %%xx
         %%include_subject_list{num_include_subject}.visit_list{1}.baseline_age = 0;
         
         visit_num_list = [visit_num_list num_include_visit];
         trace_year_list = [trace_year_list trace_year];
         total_num_include_visit = total_num_include_visit + num_include_visit;

     end

     
     
end     

include_subject_list = include_subject_list(1:num_include_subject);
total_num_include_visit
num_include_subject

str = sprintf('num subject = %d, num visits = %d\n', num_include_subject, total_num_include_visit)
str = sprintf('visit num = %f +- %f, max = %d\n', mean(visit_num_list), std(visit_num_list), max(visit_num_list))
str = sprintf('trace year = %f +- %f, min = %f, max = %f\n', mean(trace_year_list), std(trace_year_list), min(trace_year_list), max(trace_year_list))

global out_dir;

figure,
hist(visit_num_list)
title('visit count');
str = sprintf('%s/hist_visit_num.png', out_dir);
saveas(gcf, str, 'png');

figure,
hist(trace_year_list);
title('trace year');
str = sprintf('%s/trace_year.png', out_dir);
saveas(gcf, str, 'png');

for s = 1:num_include_subject    
     num_visit = include_subject_list{s}.num_visit;
     for v = 1:num_visit
        include_subject_list{s}.visit_list{v}.data = include_subject_list{s}.visit_list{v}.data(data_idx);
        
        %% add on 2015/6/2
        include_subject_list{s}.visit_time_list(v) = include_subject_list{s}.visit_list{v}.month / 12.0;
        for d = 1:3
            include_subject_list{s}.visit_data_list(v, d) = include_subject_list{s}.visit_list{v}.data(d);
        end  
     end
end
