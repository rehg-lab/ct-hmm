function [subject_list, num_visit_list, trace_year_list, data1_list, data2_list, baseline_age_list, age_list] = ...
                                                    func_extract_qual_eye_data(eye_list, is_extract_g_gs, min_num_visit, fp_log)

num_eye = size(eye_list, 1);

subject_list = cell(num_eye, 1);
num_include_subject = 0;
num_other_subject = 0;

max_num_visit = 100000;

data1_list = zeros(max_num_visit, 1);
data2_list = zeros(max_num_visit, 1);
age_list = zeros(max_num_visit, 1);

trace_year_list = zeros(num_eye, 1);
num_visit_list = zeros(num_eye, 1);
baseline_age_list = zeros(num_eye, 1);

total_num_include_visit = 0;

for e = 1:num_eye   
    
    %% check if this eye has more than 1 qualified visit
    num_visit = eye_list{e}.num_visit;    
    num_include_visit = 0;
    
    has_g_or_gs = 0;

    for v = 1:num_visit
        if (eye_list{e}.visit_list{v}.include == 1)
            num_include_visit = num_include_visit + 1;                       
            %% check if there is at least one visit has g/gs diagnosis                
            if (strcmp(eye_list{e}.visit_list{v}.Dx, 'g') == 1 || ...
                strcmp(eye_list{e}.visit_list{v}.Dx, 'gs') == 1)                
                has_g_or_gs = 1;
            end
        end
    end

    is_include_subject = 0;
    if (is_extract_g_gs == 1 && has_g_or_gs == 1)
        is_include_subject = 1;
    elseif (is_extract_g_gs == 0 && has_g_or_gs == 0)
        is_include_subject = 1;
    else
        num_other_subject = num_other_subject + 1;
    end
    
    if (num_include_visit >= min_num_visit && is_include_subject == 1)
    
        num_include_subject = num_include_subject + 1;
        s = num_include_subject;

        subject_list{s}.ID = eye_list{e}.ID;
        subject_list{s}.eye = eye_list{e}.eye;
        
        subject_list{s}.visit_list = cell(num_include_visit, 1);
        subject_list{s}.num_visit = num_include_visit;
        
        %% add on 5/10/2015
        subject_list{s}.visit_time_list = zeros(num_include_visit, 1);
        dim = 2;
        subject_list{s}.visit_data_list = zeros(num_include_visit, dim);
        
        %baseline_age_list = [baseline_age_list eye_list{e}.visit_list{1}.baseline_age];
        
        k = 0;
        for v = 1:num_visit
            if (eye_list{e}.visit_list{v}.include == 1)   % this visit has complete data          
                
                k = k+1;
                total_num_include_visit = total_num_include_visit + 1;
                
                subject_list{s}.visit_list{k} = eye_list{e}.visit_list{v};
                
                %% add on 5/10/2015
                %subject_list{s}.visit_time_list(k) = round(eye_list{e}.visit_list{v}.day / 30.0);
                subject_list{s}.visit_time_list(k) = eye_list{e}.visit_list{v}.month;
                
                subject_list{s}.visit_data_list(k,1) = eye_list{e}.visit_list{v}.data(1);
                subject_list{s}.visit_data_list(k,2) = eye_list{e}.visit_list{v}.data(2);
                                
                data1_list(total_num_include_visit) = eye_list{e}.visit_list{v}.data(1);
                data2_list(total_num_include_visit) = eye_list{e}.visit_list{v}.data(2);
                age_list(total_num_include_visit) = eye_list{e}.visit_list{v}.age;
                               
            end
        end
        
        baseline_age_list(num_include_subject) = subject_list{s}.visit_list{1}.age;                                
        num_visit_list(num_include_subject) = num_include_visit;   
        
        last_day = subject_list{s}.visit_list{num_include_visit}.day;
        first_day = subject_list{s}.visit_list{1}.day;
        
        year = (last_day - first_day)/ 365.0;  
        trace_year_list(num_include_subject) = year;
               
    end
end

%% output statistics
subject_list = subject_list(1:num_include_subject);

data1_list = data1_list(1:total_num_include_visit);
data2_list = data2_list(1:total_num_include_visit);
age_list = age_list(1:total_num_include_visit);

baseline_age_list = baseline_age_list(1:num_include_subject);
num_visit_list = num_visit_list(1:num_include_subject);
trace_year_list = trace_year_list(1:num_include_subject);

if (is_extract_g_gs == 1)
    str_included_type = 'g/gs';
else
    str_included_type = 'normal';
end 

fprintf(fp_log, 'num included subject = %d (%s)\n', num_include_subject, str_included_type);
fprintf(fp_log, 'num other subject = %d\n', num_other_subject);
fprintf(fp_log, 'num included visits = %d\n', total_num_include_visit);

