function [subject_list, trace_year_list, num_visit_list] = func_glaucoma_select_train_visit_age(age_l, age_u, ori_subject_list)

dim = 2;
global fp_log;

num_ori_subject = size(ori_subject_list, 1);
num_include_subject = 0;
subject_list = cell(num_ori_subject, 1);

trace_year_list = zeros(num_ori_subject, 1);
num_visit_list = zeros(num_ori_subject, 1);
ID_list = zeros(num_ori_subject, 1);
total_num_include_visit = 0;

for s = 1:num_ori_subject
    
    %% check if this eye has more than 1 qualified visit
    num_ori_visit = ori_subject_list{s}.num_visit;    
    num_include_visit = 0;
    include_visit_idx = [];
    
    for v = 1:num_ori_visit
        if (ori_subject_list{s}.visit_list{v}.age >= age_l && ...
            ori_subject_list{s}.visit_list{v}.age < age_u)
            num_include_visit = num_include_visit + 1;
            include_visit_idx = [include_visit_idx v];
        end
    end    
    
    if (num_include_visit >= 1) % then include this subject   
        
        num_include_subject = num_include_subject + 1;
        m = num_include_subject;
        total_num_include_visit = total_num_include_visit + num_include_visit;
                
        %% compute trace year
        first_include_idx = include_visit_idx(1);
        last_include_idx = include_visit_idx(end);
        last_month = ori_subject_list{s}.visit_list{last_include_idx}.month;
        first_month = ori_subject_list{s}.visit_list{first_include_idx}.month;
        trace_year = (last_month - first_month) / 12.0;        
        trace_year_list(m) = trace_year;   
        num_visit_list(m) = num_include_visit;

        subject_list{m}.ID = ori_subject_list{s}.ID;
        ID_list(m) = subject_list{m}.ID;
        
        subject_list{m}.eye = ori_subject_list{s}.eye;
        
        subject_list{m}.visit_list = cell(num_include_visit, 1);
        subject_list{m}.visit_time_list = zeros(num_include_visit, 1);        
        subject_list{m}.visit_data_list = zeros(num_include_visit, dim);
        
        subject_list{m}.num_visit = num_include_visit;
              
        for i = 1:num_include_visit
           k = include_visit_idx(i);           
           subject_list{m}.visit_list{i} = ori_subject_list{s}.visit_list{k};
           subject_list{m}.visit_time_list(i) = ori_subject_list{s}.visit_time_list(k);
           subject_list{m}.visit_data_list(i, :) = ori_subject_list{s}.visit_data_list(k, :);           
        end % i

    end % if
end % s

subject_list = subject_list(1:num_include_subject);
trace_year_list = trace_year_list(1:num_include_subject);
num_visit_list = num_visit_list(1:num_include_subject);
ID_list = ID_list(1:num_include_subject);

str = sprintf('num included eyes = %d, total visits = %d\n', num_include_subject, total_num_include_visit);
fprintf(fp_log, str);
fprintf(str);

str = sprintf('num of unique subject ID = %d\n', size(unique(ID_list), 1));
fprintf(fp_log, str);
fprintf(str);

str = sprintf('trace year: mean:%f, std=%f\n', mean(trace_year_list), std(trace_year_list));
fprintf(fp_log, str);
fprintf(str);

str = sprintf('num visit: mean:%f, std=%f\n', mean(num_visit_list), std(num_visit_list));
fprintf(fp_log, str);
fprintf(str);
