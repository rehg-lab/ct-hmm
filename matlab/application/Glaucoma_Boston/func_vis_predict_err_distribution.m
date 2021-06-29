function func_vis_predict_err_distribution(out_dir, overall_CTHMM_abs_err, overall_LR_abs_err, overall_bayes_LR_abs_err, overall_global_LR_abs_err)

global num_min_hist_visit;

fontsize = 6;

data_dim = size(overall_CTHMM_abs_err, 1);

%data_name = {'VFI', 'RNFL'};
global data_setting;
data_name = data_setting.type_name_ls;


for d = 1:data_dim

figure,
hist(overall_CTHMM_abs_err(d, :));
str = sprintf('Histogram of CT-HMM prediction error (%s: %.2f +- %.2f)', data_name{d}, mean(overall_CTHMM_abs_err(d, :)), std(overall_CTHMM_abs_err(d, :)));
title(str, 'FontSize', fontsize);
xlabel(data_name{d});
filename = sprintf('%s/CTHMM_pred_err_hist_%s', out_dir, data_name{d});
print(filename, '-dpng', '-r600'); %<-Save as PNG with 600 DPIs
close(gcf);

if (num_min_hist_visit >= 2)    
    figure,
    hist(overall_LR_abs_err(d, :));
    str = sprintf('Histogram of Linear Regression prediction error (%s: %.2f +- %.2f)', data_name{d}, mean(overall_LR_abs_err(d, :)), std(overall_LR_abs_err(d, :)));
    title(str, 'FontSize', fontsize);
    xlabel(data_name{d});
    filename = sprintf('%s/LR_pred_err_hist_%s', out_dir, data_name{d});
    print(filename, '-dpng', '-r600'); %<-Save as PNG with 600 DPIs
    close(gcf);

    figure,
    hist(overall_bayes_LR_abs_err(d, :));
    str = sprintf('Histogram of Bayesian Linear Regression prediction error (%s: %.2f +- %.2f)', data_name{d}, mean(overall_bayes_LR_abs_err(d, :)), std(overall_bayes_LR_abs_err(d, :)));
    title(str, 'FontSize', fontsize);
    xlabel(data_name{d});
    filename = sprintf('%s/Bayes_pred_err_hist_%s', out_dir, data_name{d});
    print(filename, '-dpng', '-r600'); %<-Save as PNG with 600 DPIs
    close(gcf);
end
    
figure,
hist(overall_global_LR_abs_err(d, :));
str = sprintf('Histogram of Global Linear Regression prediction error (%s: %.2f +- %.2f)', data_name{d}, mean(overall_global_LR_abs_err(d, :)), std(overall_global_LR_abs_err(d, :)));
title(str, 'FontSize', fontsize);
filename = sprintf('%s/globalLR_pred_err_hist_%s', out_dir, data_name{d});
print(filename, '-dpng', '-r600'); %<-Save as PNG with 600 DPIs
close(gcf);

%% draw the error difference between HMM and Bayesian, HMM and LR
if (num_min_hist_visit >= 2)
    err_diff = overall_CTHMM_abs_err(d, :) - overall_bayes_LR_abs_err(d, :);
    figure,
    hist(err_diff);
    str = sprintf('Histogram of the difference of CT-HMM and Bayesian LR prediction error (%s: %.2f +- %.2f)', data_name{d}, mean(err_diff), std(err_diff));
    title(str, 'FontSize', fontsize);
    xlabel(data_name{d});
    filename = sprintf('%s/CTHMM_BayesLR_pred_diff_hist_%s', out_dir, data_name{d});
    print(filename, '-dpng', '-r600'); %<-Save as PNG with 600 DPIs
    close(gcf);

    figure,
    err_diff = overall_CTHMM_abs_err(d, :) - overall_LR_abs_err(d, :);
    hist(err_diff);
    str = sprintf('Histogram of the difference of CT-HMM and LR prediction error (%s: %.2f +- %.2f)', data_name{d}, mean(err_diff), std(err_diff));
    title(str, 'FontSize', fontsize);
    xlabel(data_name{d});
    filename = sprintf('%s/CTHMM_LR_pred_diff_hist_%s', out_dir, data_name{d});
    print(filename, '-dpng', '-r600'); %<-Save as PNG with 600 DPIs
    close(gcf);
end

end
