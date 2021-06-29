%% decode_simulation test
%% test decoding accuracy under different noise levels

out_directory = 'test_output';  %% set up output directory

addpath('decode/baseline_SSA');
addpath('decode');
addpath('precompute');
addpath('likelihood');
addpath('common');
addpath('common/fast_expm_A_t');
addpath('MD/MD_model/');

%% global variables
global Q_mat;
global syn_Q_mat;

global out_dir;

%% observation sequences 
global obs_seq_list; % contain the generated synthetic data

%% setup output
mkdir(out_directory);

%% log file handle
global fp_log;
str = sprintf('%s/decode_log.txt', out_directory);
fp_log = fopen(str, 'wt');

%% setup synthetic data
syn_qi_range = [1 5]; % set qi to be in range [1 5], except the absorbing state
run_list = [1]'; % 1 random runs
num_run = length(run_list);

%% set up the testing state space
num_dim = 1;
num_state_per_dim = 5; % a 5-state model
neighbor_setting = 4; % 4 for fully-connected model(1: forward link only, 2: backward link only, 3: both directions, 4: fully connected)

%% synthetic data setting
syn_data_config_list = [1];
num_syn_config = length(syn_data_config_list);
test_method_list = [1];
num_test_method = length(test_method_list);
num_decoding_accu_result = 2; % one for continuous-time result, one for observation-time result

all_result = zeros(num_run, num_syn_config, num_test_method, num_decoding_accu_result);

%% for each random run
for r = 1:num_run
        
    %% specify the seed for random number
    run_idx = run_list(r);
    sd = run_idx * 100;
    rng(sd);
    
    str = sprintf('===== Run Idx %d: (seed = %d) ======\n', run_idx, sd);
    CTHMM_print_log(str);
    
    %% create state list
    state_sigma = 0.25; %state sigma is set to be 0.25 of the state's data range
    CTHMM_sim_create_MD_syn_state_list(num_dim, num_state_per_dim, state_sigma);
     
    %% create synthetic Q mat and init learning Q mat    
    CTHMM_sim_gen_syn_Q_mat(syn_qi_range, neighbor_setting);
    Q_mat = syn_Q_mat;
    
    %% create state reachability mat
    CTHMM_precompute_state_reach_mat();
                 
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
        
        %% generate observation sequences
        CTHMM_sim_batch_gen_syn_obs_seq(syn_data_settings); %% need to record the underlying true state sequence for the decoding experiment
        
        %% start the decoding
        num_obs_seq = size(obs_seq_list, 1);        
          
        %% for each testing method
        for m = 1:num_test_method
            
            method = test_method_list(m);
            method_name = method_name_list{method};
            
            str = sprintf('===== Method %d: %s =====\n', method, method_name);
            CTHMM_print_log(str);                                 
            
            %% do decoding
            has_ground_truth = 1;
            [ave_conti_state_decode_err, ave_obstime_state_decode_err] = CTHMM_decode_batch(obs_seq_list, expect_tau_method, has_ground_truth);            
            
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
        all_obstime_state_decode_err = all_result(:, s, m, 1);
        
        str = sprintf('conti_state_decode_err = %f +- %f\n', mean(all_conti_state_decode_err), std(all_conti_state_decode_err));
        CTHMM_print_log(str); 
        
        str = sprintf('obstime_state_decode_err = %f +- %f\n', mean(all_obstime_state_decode_err), std(all_obstime_state_decode_err));
        CTHMM_print_log(str);                
    end
end

fclose(fp_log);
