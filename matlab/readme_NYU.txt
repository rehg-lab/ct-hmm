* Datasheet format (use .xlsx file)
  
  For RNFL uncalibrated datasheet (such as "longitudinal dataset new.xlsx" where we used in NIPS), the RNFL data is the raw data from different OCT machines, my program does explicit RNFL data calibration). The program will parse the following fields:   
  'subject', 'eye', 'daysfrombaseline', 'Dx', 'baseline_age', 'age', 'vf', 'VFI', 'rnfl'	

  For RNFL calibrated datasheet (the datasheets that NYU used recently),  
  my program will parse the following fields:  
  'subject', 'eye', 'daysfrombaseline', 'Dx', 'baseline_age', 'age', 'data1', 'data2'
  where 'data1' is for your dimention 1, and 'data2' for your dimenion 2.
  (i.e., if you want RNFL data as your dimension 1, you might copy RNFL data to another column and with field name "data1"). No explicit calibration is needed in the program.
  
* For the following scripts explanation, the scripts with "calibrated" suffix will load calibrated datasheet, while the versions without this string are the old codes we used for the uncalibrated data (NIPS paper).

--------------
To train glaucoma models, go to application\Glaucoma_Boston:
  
* To train and visualize the learned CT-HMM models, run 
  "run_glaucoma_age_group_visualization_calibrated" in Matlab command line:
  
  - Put your datasheet under Glaucoma_Boston_parse\  
  - Change to your datasheet name
  - Output is in output_vis_minvisit5\ 
    -- visualization of state transition chart is in age_0_120\CV_1\IterXX\vis_Nij_mat.png  (XX is the last iteration) 
    -- State count and learning time information is in age_0_120\CV_1\log.txt
	-- Data statistics is in log.txt and data_stat\
	-- Note: during learning, if using the hybrid method (learn_method=4), the first iteration will be slower (about 15 minutes, as EXPM method is used), but the remaining iterations will be much faster (each <1 min, used eigen method)	   
  - You may change your state definition
  - You may learn models for different age groups
  
* To predict future observations for each eye (use 10-fold cross validation), run
  "run_glaucoma_prediction_cross_validation_calibrated"
  
  - Change to your datasheet name
  - In testing, for each testing eye's data, the first 4 visits were used as history, and the system will predict the values for the remaining future data (you may modify the variable num_min_hist_visit)
  - Output is in output_predict_cv10_minvisit5\ 
    Prediction result is in the end of log_10run.txt
  - The total running time for 10-fold cross validation may take hours	
  
* To decode and visualize the underlying continuous state path for each eye's longitudinal data, run
  "run_glaucoma_continuous_decoding_calibrated"
  
  - Change to your datasheet name
  - The default in this function will train a model and do decoding for each included eye
  - If you have trained a model, you may specify a pretrained CT-HMM model (without training again)
  - Output is in output_decoding\ 
  - Note: the decoding algorithm has not been formally published yet. 
    We have a paper draft for the method 
  
  For decoding, there are two options in my code.
  Option 1: you may train the model again using all the data, and then do the decoding. To do so, set run_model_training = 1
  Option 2: if you already train the model (by using run_glaucoma_age_group_visualization_calibrated.m), then you can just set the path for the trained model. For example, you may set top_out_dir = 'output_vis_minvisit5/age_0_120'

-------------
Reference paper: 
"Efficient Learning of Continuous-Time Hidden Markov Models for Disease Progression", NIPS 2015
https://papers.nips.cc/paper/5792-efficient-learning-of-continuous-time-hidden-markov-models-for-disease-progression
 