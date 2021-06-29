%% do complete path decoding using ViterbiExpect method
function [overall_state_seq, overall_state_dur_seq, outer_best_state_seq] = CTHMM_decode_complete_path_viterbiexpect(obs_seq, expect_dur_method)

max_state_seq_len = 100; %% can set to higher when the length exceeds
overall_state_seq = zeros(max_state_seq_len, 1);
overall_state_dur_seq = zeros(max_state_seq_len, 1);
cur_state_seq_len = 0;

eigen_num_fail = 0;
%closed_num_fail = 0; %% closed-form method

%% outer state decoding using Viterbi
[outer_best_state_seq, outer_dur_seq, best_log_prob, log_Pt_list] = CTHMM_decode_outer_viterbi(obs_seq);

num_visit = obs_seq.num_visit;

global cached_inner_best_state_seq;
global is_decode_simulation;
global is_glaucoma_decoding;


for v = 1:(num_visit-1)
    
    %v
    
    %% two interval between two successive observations
    %T = outer_dur_seq(v+1) - outer_dur_seq(v);
    T = outer_dur_seq(v);
    
    %% For every two successive outer decoded states, do inner state decoding using the method from "state sequence analysis" paper (molecular)    
    start_s = outer_best_state_seq(v);
    end_s = outer_best_state_seq(v+1);
    
    if (is_decode_simulation == 0)
        
        v
        
        %% compute best inner state sequence
        
        if (is_glaucoma_decoding == 1 && (start_s == end_s))
            inner_best_state_seq = start_s;
            inner_best_state_seq
        else
            tic;
            [inner_best_state_seq] = CTMC_decode_most_probable_state_seq_SSA(start_s, end_s, T);
            tEnd = toc;
            str = sprintf('SSA: %d min, %f sec\n', floor(tEnd/60),rem(tEnd,60))

            inner_best_state_seq
        end
        
        
        
        %% compute best inner duration sequence
        %% For the computed inner state path, compute the expected state dwell time for each inner state  
        tic;
        
        if (expect_dur_method == 1)    % Expm

            
            %% expm method
            [dur_list] = CTHMM_decode_expected_dur_for_a_path_Expm(inner_best_state_seq, T);    
            
            
            
        elseif (expect_dur_method == 2)   % Unif 

            %% unif method        
            [dur_list] = CTHMM_decode_expected_dur_for_a_path_Unif(inner_best_state_seq, T);

        elseif (expect_dur_method == 3) % Eigen

            %% eigen method        
            [dur_list, is_success] = CTHMM_decode_expected_dur_for_a_path_Eigen(inner_best_state_seq, T);        
            if (is_success == 0)
                eigen_num_fail = eigen_num_fail + 1;       
            end

        elseif (expect_dur_method == 4)
    %         %% closed-form method        
    %         [dur_list, is_success, a_table, b_table, c_table] = CTHMM_decode_expected_dur_for_a_path_closedform(inner_best_state_seq, T);                
    %         if (is_success == 0)
    %             closed_num_fail = closed_num_fail + 1;
    %         end
        else       
            disp('unknwn expect_dur_method');
        end
        
        tEnd = toc;
        str = sprintf('Compute Dur: %d min, %f sec\n', floor(tEnd/60),rem(tEnd,60))
         
    else % is decoding simulation, use precomputation and caching
            
        %% just for speeding up the simulation experiment
        if (isempty(cached_inner_best_state_seq{start_s, end_s}.state_ls)==0) 
            inner_best_state_seq = cached_inner_best_state_seq{start_s, end_s}.state_ls; % here, we assume that T is always the same (in the simulation test)        
            dur_list = cached_inner_best_state_seq{start_s, end_s}.dur_ls;
        else
            disp('please precompute the best state sequence');
            exit;
        end        
    end
            
    %% if needed, enlarge the size of the overall state sequence
    num_inner_state = length(inner_best_state_seq);        
    if ((cur_state_seq_len + length(inner_best_state_seq)) > max_state_seq_len) %% enlarge the overall state sequence size
        temp = overall_state_seq;
        max_state_seq_len = max_state_seq_len * 2;
        overall_state_seq = zeros(max_state_seq_len, 1);
        overall_state_seq(1:cur_state_seq_len) = temp(1:cur_state_seq_len);
    end
        
    %% append the inner_best_state_seq to the overall_state_seq, and dur_list to the overall_state_dur_seq
    if (cur_state_seq_len > 0 && (inner_best_state_seq(1) == overall_state_seq(cur_state_seq_len)))
        overall_state_seq(cur_state_seq_len:(cur_state_seq_len+num_inner_state-1)) = inner_best_state_seq;
        overall_state_dur_seq(cur_state_seq_len) = overall_state_dur_seq(cur_state_seq_len) + dur_list(1);
        if (num_inner_state > 1)            
            overall_state_dur_seq((cur_state_seq_len+1):(cur_state_seq_len+num_inner_state-1)) = dur_list(2:end);
        end
        cur_state_seq_len = cur_state_seq_len + num_inner_state - 1;
    else %% direct add the new sequence
        overall_state_seq((cur_state_seq_len+1):(cur_state_seq_len+num_inner_state)) = inner_best_state_seq;
        overall_state_dur_seq((cur_state_seq_len+1):(cur_state_seq_len+num_inner_state)) = dur_list;
        cur_state_seq_len = cur_state_seq_len + num_inner_state;        
    end          
    
end

overall_state_seq = overall_state_seq(1:cur_state_seq_len);
overall_state_dur_seq = overall_state_dur_seq(1:cur_state_seq_len);
