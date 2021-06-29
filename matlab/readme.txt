Simulation test:
- run CTHMM_test_script.m to do a 5-state simulation test for the specified learning methods 

Before running the script, please change your current Matlab path to the path where you put the CT-HMM project. You can use cd('your CT-HMM path') command in Matlab command window.

For the simulation, data will be generated and put in the global variable obs_seq_list. If you want to load your own data, fill the data fields in obs_seq_list before learning:
obs_seq_list{subject_idx}.num_visit
obs_seq_list{subject_idx}.visit_time_list(visit_idx)
obs_seq_list{subject_idx}.visit_data_list(visit_idx, data_dim)


Note: 
* we have recently developed a faster eigen method (one order of S faster than algo3 in the NIPS supplemental, and it becomes the fastest algorithm among all the three alternatives. This will be published in a book chapter in August 2017.
* If Eigen method failed in the first iteration with uniform initialization, you may try to add a small random noise to the Q_mat initialization and try again
* If Eigen failed in any iteration, in practice one can also switch to expm/unif method for that iteration and then switch back to eigen method (you need to modify the code for switching by yourself for this feature.)

