global state_list;
state_list = cell(36, 1);

global is_use_individual_Q_mat;
is_use_individual_Q_mat = 0;

global is_use_distinct_time_grouping;
is_use_distinct_time_grouping = 0;

global Q_mat;
% Q_mat = rand(2,2);
Q_mat = load('test10Q.txt');

global state_init_prob_list;
state_init_prob_list = ones(100)./100;

obs_seq.num_visit = 5;
obs_seq.visit_time_list = [25.742465753424696,29.227397260273968,33.43561643835616,38.005479452054786,41.45753424657535];
obs_seq.visit_data_list = [[0.0, 0.0], [0.0, 0.0], [9.0, 5.0], [32.0, 18.0], [34.0, 4.0]];
obs_seq.has_compute_data_emiss_prob = 1;

% time visit, state
obs_seq.log_data_emiss_prob_list = load("test10_logemission.txt");

[outer_best_state_seq, outer_dur_seq, best_log_prob, log_Pt_list] = CTHMM_decode_outer_viterbi(obs_seq);

disp(outer_best_state_seq)
disp(outer_dur_seq)

for v = 2
    
    %v
    
    %% two interval between two successive observations
    %T = outer_dur_seq(v+1) - outer_dur_seq(v);
    T = outer_dur_seq(v);
    
    %% For every two successive outer decoded states, do inner state decoding using the method from "state sequence analysis" paper (molecular)    
    start_s = outer_best_state_seq(v);
    end_s = outer_best_state_seq(v+1);
    
    [inner_best_state_seq, best_prob_SSA] = CTMC_decode_most_probable_state_seq_SSA(start_s, end_s, T);
    disp("inner state sequence")
    disp(start_s)
    disp(end_s)
    disp(inner_best_state_seq)
    disp(best_prob_SSA)
end
