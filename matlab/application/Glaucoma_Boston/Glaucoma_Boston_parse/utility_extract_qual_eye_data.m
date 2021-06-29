function [subject_list, num_visit_list, trace_year_list, VFI_list, RNFL_list, rimarea_list, baseline_age_list, age_list] = ...
utility_extract_qual_eye_data(eye_list, is_extract_g_gs, data_idx_include_list,  fp_log)

num_eye = size(eye_list, 1);

subject_list = cell(num_eye, 1);
num_include_subject = 0;
num_other_subject = 0;

VFI_list = [];
RNFL_list = [];
rimarea_list = [];

total_num_include_visit = 0;
trace_year_list = [];
num_visit_list = [];
baseline_age_list = [];
age_list = [];

for e = 1:num_eye   
    
    %% check if this eye has more than 1 qualified visit
    num_visit = eye_list{e}.num_visit;    
    num_include_visit = 0;
    
    has_g_or_gs = 0;

    for v = 1:num_visit
        
        
        eye_list{e}.visit_list{v}.include = 0;
        
        if (data_idx_include_list(1) == 1) % include VFI
            if (eye_list{e}.visit_list{v}.vf == 0 || isnan(eye_list{e}.visit_list{v}.VFI) == 1)
                continue;
            end            
        end
        
        if (data_idx_include_list(2) == 1) % include RNFL
            if (eye_list{e}.visit_list{v}.has_qual_OCT == 0) %isnan(eye_list{e}.visit_list{v}.norm_rnfl_stratus) == 1)
                continue;
            end            
        end
        
        if (data_idx_include_list(3) == 1) % include rim
            if (isnan(eye_list{e}.visit_list{v}.rimarea) == 1)
                continue;
            end            
        end
        
        eye_list{e}.visit_list{v}.include = 1;
        num_include_visit = num_include_visit + 1;                       
            
        %% check if there is at least one visit has g/gs diagnosis                
        if (strcmp(eye_list{e}.visit_list{v}.Dx, 'g') == 1 || ...
            strcmp(eye_list{e}.visit_list{v}.Dx, 'gs') == 1)                
            has_g_or_gs = 1;
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
    
    if (num_include_visit >= 1 && is_include_subject == 1)
    
        num_include_subject = num_include_subject + 1;
        s = num_include_subject;

        subject_list{s}.ID = eye_list{e}.ID;
        subject_list{s}.eye = eye_list{e}.eye;
        
        subject_list{s}.visit_list = cell(num_include_visit, 1);
        subject_list{s}.num_visit = num_include_visit;
        baseline_age_list = [baseline_age_list eye_list{e}.visit_list{1}.baseline_age];
        
        k = 0;
        for v = 1:num_visit
            
            if (eye_list{e}.visit_list{v}.include == 1)                
            
                k = k+1;
                
                subject_list{s}.visit_list{k} = eye_list{e}.visit_list{v};  
                
                subject_list{s}.visit_list{k}.data = [];
                for d = 1:3 %%%%
                    if (data_idx_include_list(d) == 1)
                        subject_list{s}.visit_list{k}.data = [subject_list{s}.visit_list{k}.data eye_list{e}.visit_list{v}.data(d)];
                    end                    
                end
                
                
                VFI_list = [VFI_list eye_list{e}.visit_list{v}.data(1)];
                RNFL_list = [RNFL_list eye_list{e}.visit_list{v}.data(2)];                
                rimarea_list = [rimarea_list eye_list{e}.visit_list{v}.data(3)];
                age_list = [age_list eye_list{e}.visit_list{v}.age];
                                
            end
            
        end
        
        last_day = subject_list{s}.visit_list{num_include_visit}.day;
        first_day = subject_list{s}.visit_list{1}.day;
        
        duration_day = (last_day - first_day) + 1;        
        year = duration_day / 365.0;
        
        num_visit_list = [num_visit_list num_include_visit];
        trace_year_list = [trace_year_list year];
        total_num_include_visit = total_num_include_visit + num_include_visit;
        
    end
end

%% output statistics
subject_list = subject_list(1:num_include_subject);

num_include_subject
num_other_subject
total_num_include_visit

fprintf(fp_log, 'num included subject = %d\n', num_include_subject);
fprintf(fp_log, 'num other subject = %d\n', num_other_subject);
fprintf(fp_log, 'num include visits = %d\n', total_num_include_visit);

