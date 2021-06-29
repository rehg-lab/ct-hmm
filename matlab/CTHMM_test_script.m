%% 5-state simulation experiments

%% set up simulation
out_directory = '../simulation_5state_syn1_m1';
syn_qi_range = [1 5]; % set qi to be in range [1 5], except the absorbing state
run_list = [1:1:5]'; % do 5 random runs
%run_list = [1]'; % 1 random runs

syn_data_config_list = [1]'; % 10^5 visits, sampling rate: smallest holding time * 0.5; write the config in CTHMM_sim_set_syn_data_config(syn_config_idx)
test_method_list = [1]'; %test 1:soft(expm), 3:soft(unif), 5:soft(new Eigen, one order of S faster than the algo we published in NIPS 2005),  2:hard(expm), 4:hard(unif), 6:hard(Eigen)

%test_method_list = [1:1:6];

%% Note: we have developed a faster eigen method (one order of S faster than algo3 in the NIPS supplemental, and it becomes the fastest algorithm among all 3 alternatives
%% If Eigen method failed in the first iteration, you may try to add a small random noise to the Q_mat initialization and try again
%% If Eigen failed in any iteration, in practice one can switch to expm/unif method for that iteration and then switch back to eigen method
                         
%% set up state space
num_dim = 2;
num_state_per_dim = 3; % a 5-state model
neighbor_setting = 4; % 4: fully-connected model (means that each state connects to each other state) (1: forward link only, 2: backward link only, 3: both foward and backward directions)

%% run
CTHMM_learn_simulation(out_directory, syn_qi_range, test_method_list, syn_data_config_list, run_list, num_dim, num_state_per_dim, neighbor_setting);
