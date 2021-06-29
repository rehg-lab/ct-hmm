all_variable = load('output_predict_cv10_minvisit5_PittNew_tol10-6_def7/result_variables');

addpath('../../prediction');

global fp_log;
fp_log = fopen('redo_state.txt', 'wt');


overall_CTHMM_abs_err = all_variable.overall_CTHMM_abs_err;
overall_LR_abs_err = all_variable.overall_LR_abs_err;
overall_global_LR_abs_err = all_variable.overall_global_LR_abs_err;
overall_bayes_LR_abs_err = all_variable.overall_bayes_LR_abs_err;
func_pred_result_statistic_test(overall_CTHMM_abs_err, overall_LR_abs_err, overall_global_LR_abs_err, overall_bayes_LR_abs_err);

fclose(fp_log);