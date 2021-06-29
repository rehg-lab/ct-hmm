function [overall_conti_state_decode_err, overall_obs_time_state_decode_err] = CTHMM_decode_batch(obs_seq_list, expect_tau_method, has_ground_truth)

global is_use_distinct_time_grouping;
is_use_distinct_time_grouping = 1;

global state_list;

if (is_use_distinct_time_grouping == 1)
    CTHMM_precompute_distinct_time_Pt_list();
end
    
num_obs_seq = size(obs_seq_list, 1);

conti_state_decode_err_list = zeros(num_obs_seq, 1);
obs_state_decode_err_list = zeros(num_obs_seq, 1);

accum_diff_time = 0.0;
accum_diff_outer_state = 0;
accum_obs_time = 0.0;
accum_num_visit = 0;

accum_direct_diff_outer_state = 0;

for i = 1:num_obs_seq

    i
    
    %% the main decoding function
    obs_seq_list{i}.has_compute_data_emiss_prob = 0;
    [decoded_conti_state_seq, decoded_conti_dwelltime_seq, decoded_obs_state_seq] = CTHMM_decode_complete_path_viterbiexpect(obs_seq_list{i}, expect_tau_method);
    
    len_decoded_conti_state_seq = length(decoded_conti_state_seq);
    str = sprintf('len_decoded_conti_state_seq = %d\n', len_decoded_conti_state_seq)
    
    
    num_visit = obs_seq_list{i}.num_visit;
    total_obs_time = obs_seq_list{i}.visit_time_list(num_visit) - obs_seq_list{i}.visit_time_list(1);
    
    accum_num_visit = accum_num_visit + num_visit;
    accum_obs_time = accum_obs_time + total_obs_time;
    
    %% compute accuracy
    if (has_ground_truth == 1)                
        
        %% compute continuous-time decoding accuracy
        num_true_conti_state = length(obs_seq_list{i}.ori_state_chain.state_idx_list);        
        
        total_diff_time = 0.0; % how much time the two state sequences are different
        state_begin_time = 0.0;
        
        % for each state in the ground truth
        for t = 1:num_true_conti_state
            
            true_state = obs_seq_list{i}.ori_state_chain.state_idx_list(t);                       
            state_end_time = state_begin_time + obs_seq_list{i}.ori_state_chain.state_dur_list(t);
            if (state_end_time > total_obs_time)
                state_end_time = total_obs_time;
            end
            
            %% query the state sequence of this time segment (begin_time ~ end_time)
            [seg_decoded_state_seq, seg_decoded_dwelltime_seq] = CTHMM_decode_query_state_seq_by_time(decoded_conti_state_seq, decoded_conti_dwelltime_seq, state_begin_time, state_end_time);
            
            %% here, just check that sum of seg_decoded_dwelltime_seq equals the total time interval (state_end_time - state_begin_time)
            sum_seg_time = sum(seg_decoded_dwelltime_seq);
            query_time = state_end_time - state_begin_time;
            if (abs(sum_seg_time - query_time) > 0.000000001)
                str = sprintf('ERR in function CTHMM_decode_query_state_seq_by_time')
                CTHMM_print_log(str);
            end
            
            num_decoded_state = length(seg_decoded_state_seq);
            
            %% compute decoding accuracy (in continuous time)
            for d = 1:num_decoded_state                
                decoded_state = seg_decoded_state_seq(d);
                if (decoded_state ~= true_state)
                    total_diff_time = total_diff_time + seg_decoded_dwelltime_seq(d);
                end
            end
            
            %% set new begin time
            state_begin_time = state_end_time;
            if (state_begin_time >= total_obs_time)
                break;
            end
        end
        
        %% compute overall continuous time state decoding accuracy        
        conti_state_decode_err_list(i) = total_diff_time / total_obs_time;
        %% accum for all sequences
        accum_diff_time = accum_diff_time + total_diff_time;        
        
        %% now compute the decoding accuracy at the observation time (outer viterbi decoding)           
        %% obs_seq_list{seq}.visit_true_state_list(v): this record the true state at observation time for this chain        
        num_diff_outer_state = 0;
        for v = 1:num_visit
            if (obs_seq_list{i}.visit_true_state_list(v) ~= decoded_obs_state_seq(v))
                num_diff_outer_state = num_diff_outer_state + 1;
            end
        end
        obs_state_decode_err_list(i) = double(num_diff_outer_state) / double(num_visit);
        
        %% accum for all sequences
        accum_diff_outer_state = accum_diff_outer_state + num_diff_outer_state;
        
        %% direct observation mapping
        for v = 1:num_visit
            obs_data = obs_seq_list{i}.visit_data_list(v);
            true_state_idx = obs_seq_list{i}.visit_true_state_list(v);
            state_mu = state_list{true_state_idx}.mu;            
            if (abs(state_mu - obs_data) > 0.5)
                accum_direct_diff_outer_state = accum_direct_diff_outer_state + 1;
            end
        end        
        
    end
    
end

seq_ave_conti_state_decode_err = mean(conti_state_decode_err_list)
seq_ave_obs_time_state_decode_err = mean(obs_state_decode_err_list)

overall_conti_state_decode_err = accum_diff_time / accum_obs_time
overall_obs_time_state_decode_err =  accum_diff_outer_state / accum_num_visit

overall_direct_diff_outer_state_rate = accum_direct_diff_outer_state / accum_num_visit