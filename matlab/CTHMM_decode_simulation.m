%% decode_simulation test
%% test decoding accuracy under different noise levels

global out_dir;
out_dir = 'test_output';  %% set up output directory

addpath('decode/baseline_SSA');
addpath('decode');
addpath('decode/CTMC_decode_most_probable_state_seq_SSA.m');
addpath('precompute');
addpath('likelihood');
addpath('common');
addpath('common/fast_expm_A_t');
addpath('MD/MD_model/');
addpath('simulation');

%% global variables
global Q_mat;
global syn_Q_mat;

%% observation sequences 
global obs_seq_list; % contain the generated synthetic data

%% setup output
mkdir(out_dir);

%% log file handle
global fp_log;
str = sprintf('%s/decode_log.txt', out_dir);
fp_log = fopen(str, 'wt');

%% setup synthetic data
syn_qi_range = [1 5]; % set qi to be in range [1 5], except the absorbing state
%run_list = [1:1:5]'; % 5 random runs
run_list = [1]'; % 1 random runs

num_run = length(run_list);

%% set up the testing state space
num_dim = 1;
num_state_per_dim = 5; % a 5-state model
neighbor_setting = 4; % 4 for fully-connected model(1: forward link only, 2: backward link only, 3: both directions, 4: fully connected)

%% synthetic data setting
%syn_data_config_list = [20:27];
syn_data_config_list = [20:23];

num_syn_config = length(syn_data_config_list);
global syn_data_settings;
global state_sigma;

method_name_list = {'Expm', 'Unif', 'Eigen'};
test_method_list = [1];
num_test_method = length(test_method_list);

num_decoding_accu_result = 2; % one for continuous-time result, one for observation-time result

all_result = zeros(num_run, num_syn_config, num_test_method, num_decoding_accu_result);

global is_use_individual_Q_mat;
is_use_individual_Q_mat = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% implement the cache function to speed up the simulation

global state_list;
num_state = length(state_list);

global is_decode_simulation;
global cached_inner_best_state_seq;
is_decode_simulation = 1;

global is_use_distinct_time_grouping;    
global distinct_time_list;
global learn_method;
global is_outer_soft;
is_outer_soft = 0; % user Viterbi decoding for the outer observation sequences
learn_method = 1; % expm
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% for each random run
for r = 1:num_run
        
    %% specify the seed for random number
    run_idx = run_list(r);
    sd = run_idx * 100;
    rng(sd);
    
    str = sprintf('===== Run Idx %d: (seed = %d) ======\n', run_idx, sd);
    CTHMM_print_log(str);
    
    %% for each synthetic data config
    for s = 1:num_syn_config
        
        syn_config_idx = syn_data_config_list(s);        
        str = sprintf('===== Syn Config Idx %d: =====\n', syn_config_idx);
        CTHMM_print_log(str);
                        
        %% create output dir for this setting
        %str = sprintf('%s/r%d_syn%d', out_directory, run_idx, syn_config_idx);
        %out_dir = str;
        %mkdir(out_dir);
        
        %% set up synthetic config
        CTHMM_sim_set_syn_data_config(syn_config_idx);
        
        %% create state list
        %state_sigma = 0.25; %state sigma is set to be 0.25 of the state's data range
        CTHMM_sim_create_MD_syn_state_list(num_dim, num_state_per_dim, state_sigma);

        %% create synthetic Q mat and init learning Q mat    
        CTHMM_sim_gen_syn_Q_mat(syn_qi_range, neighbor_setting);
        Q_mat = syn_Q_mat;

        %% create state reachability mat
        CTHMM_precompute_state_reach_mat();        
                
        %% generate observation sequences
        CTHMM_sim_batch_gen_syn_obs_seq(syn_data_settings); %% need to record the underlying true state sequence for the decoding experiment
        
        %% start the decoding
        num_obs_seq = size(obs_seq_list, 1);
               
        %=====================================================================================================                
        %% compute the inner state sequence between every two possible end-states, assume only one time interval
        if (is_decode_simulation == 1) %% just for speeding up the simulation experiment         
            str = sprintf('Precomputing inner best state sequences...')

            cached_inner_best_state_seq = cell(num_state, num_state);
            T = obs_seq_list{1}.visit_time_list(2) - obs_seq_list{1}.visit_time_list(1);            

            if (is_use_distinct_time_grouping == 1)        
                distinct_time_list = T;
                CTHMM_precompute_distinct_time_Pt_list();
            end

            
            SSA_runtime_list = zeros(num_state * num_state, 1);
            Best_SSA_prob_list = zeros(num_state * num_state, 1);
            Expm_runtime_list = zeros(num_state * num_state, 1);
            
            idx = 0;
            
            %% precompute the inner best state sequences
            for start_s = 1:num_state
                for end_s = 1:num_state 

                    
                    idx = idx + 1;
                    
                    str = sprintf('start_s = %d, end_s = %d\n', start_s, end_s)
            
                    %% compute best inner state sequence
                    
                    tic;
                    [inner_best_state_seq, best_prob_SSA] = CTMC_decode_most_probable_state_seq_SSA(start_s, end_s, T); 
                    tEnd = toc;
                    str = sprintf('SSA: %d min, %f sec\n', floor(tEnd/60),rem(tEnd,60));
                    inner_best_state_seq
                    best_prob_SSA
                    CTHMM_print_log(str);
                    SSA_runtime_list(idx) = tEnd;
                    Best_SSA_prob_list(idx) = best_prob_SSA;
                    
                    tic;
                    [dur_list] = CTHMM_decode_expected_dur_for_a_path_Expm(inner_best_state_seq, T);  
                    tEnd = toc;
                    str = sprintf('Expm_dur: %d min, %f sec\n', floor(tEnd/60),rem(tEnd,60));
                    CTHMM_print_log(str);
                    Expm_runtime_list(idx) = tEnd;
                    
                    %% check whether dur_list sum to T
                    sum_dur = sum(dur_list);
                    if (abs(sum_dur - T) > 0.000000001)
                        str = sprintf('ERR: dur_list not sum to T')             
                        [seg_decoded_state_seq, seg_decoded_dwelltime_seq] = CTHMM_decode_query_state_seq_by_time(decoded_conti_state_seq, decoded_conti_dwelltime_seq, state_begin_time, state_end_time);
                    end
                    
                    cached_inner_best_state_seq{start_s, end_s}.state_ls = inner_best_state_seq;
                    cached_inner_best_state_seq{start_s, end_s}.dur_ls = dur_list;
                end
            end  
            
            %% compute the average runtime and average best probability
            mean_SSA_runtime = mean(SSA_runtime_list);
            str = sprintf('AVE SSA run-time: %d min, %f sec\n', floor(mean_SSA_runtime/60),rem(mean_SSA_runtime,60))
            CTHMM_print_log(str);
            
            mean_Expm_runtime = mean(Expm_runtime_list);
            str = sprintf('AVE Expm-dur run-time: %d min, %f sec\n', floor(mean_Expm_runtime/60),rem(mean_Expm_runtime,60))
            CTHMM_print_log(str);
            
            mean_best_SSA_prob = mean(Best_SSA_prob_list);
            std_best_SSA_prob = std(Best_SSA_prob_list);
            str = sprintf('AVE best SSA prob: %f +- %f\n', mean_best_SSA_prob, std_best_SSA_prob)
            CTHMM_print_log(str);
            
        end        
        %=====================================================================================================
              
        T
        
        %% for each testing method (the method for computing the expected duration of each state)
        for m = 1:num_test_method
            
            method = test_method_list(m);
            method_name = method_name_list{method};
            
            str = sprintf('===== Method %d: %s =====\n', method, method_name);
            CTHMM_print_log(str);                                 
            
            %% do decoding
            has_ground_truth = 1;
            
            tic;
            [ave_conti_state_decode_err, ave_obstime_state_decode_err] = CTHMM_decode_batch(obs_seq_list, method, has_ground_truth);            
            tEnd = toc;
            str = sprintf('CTHMM_decode_batch: %d min, %f sec\n', floor(tEnd/60),rem(tEnd,60));
            CTHMM_print_log(str);
                    
            %% write the result
            str = sprintf('ave_conti_state_decode_err = %f, ave_obstime_state_decode_err = %f\n\n', ave_conti_state_decode_err, ave_obstime_state_decode_err);
            CTHMM_print_log(str);
            
            all_result(r, s, m, 1) = ave_conti_state_decode_err;
            all_result(r, s, m, 2) = ave_obstime_state_decode_err;
            
            %% save all_result so far
            str = sprintf('%s/all_result', out_dir);   
            save(str, 'all_result');

        end % method
                
    end % syn_config
end % run

%% save all decoding result
str = sprintf('%s/all_result', out_dir);   
save(str, 'all_result');

%CTHMM_sim_report_learn_result(all_result, test_method_list, syn_data_config_list, run_list);

%% calculate average results of all runs and print to log
for s = 1:num_syn_config
        
    syn_config_idx = syn_data_config_list(s);        
    str = sprintf('===== Syn Config Idx %d: =====\n', syn_config_idx);
    CTHMM_print_log(str);

    for m = 1:num_test_method
            
        method = test_method_list(m);
        method_name = method_name_list{method};

        str = sprintf('===== Method %d: %s =====\n', method, method_name);
        CTHMM_print_log(str);                                 
        
        all_conti_state_decode_err = all_result(:, s, m, 1);
        all_obstime_state_decode_err = all_result(:, s, m, 2);
        
        str = sprintf('conti_state_decode_err = %f +- %f\n', mean(all_conti_state_decode_err), std(all_conti_state_decode_err));
        CTHMM_print_log(str); 
        
        str = sprintf('obstime_state_decode_err = %f +- %f\n', mean(all_obstime_state_decode_err), std(all_obstime_state_decode_err));
        CTHMM_print_log(str);                
    end
end

fclose(fp_log);
