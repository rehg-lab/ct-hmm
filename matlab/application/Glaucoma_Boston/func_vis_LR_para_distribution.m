function func_vis_glaucoma_LR_para_distribution(out_dir, seq_hist_LR_para, seq_hist_bayesian_LR_para, seq_global_LR_para)

%% draw the distribution of the slope and intercept of the two methods
LR_method_name = 'Linear Regression';
LR_method_filename = 'LR';
vis_LR_para_distribution(out_dir, seq_hist_LR_para, LR_method_name, LR_method_filename);
LR_method_name = 'Bayesian Linear Regression';
LR_method_filename = 'BayesLR';
vis_LR_para_distribution(out_dir, seq_hist_bayesian_LR_para, LR_method_name, LR_method_filename);
LR_method_name = 'Global Linear Regression';
LR_method_filename = 'GlobalLR';
vis_LR_para_distribution(out_dir, seq_global_LR_para, LR_method_name, LR_method_filename);


    function vis_LR_para_distribution(out_dir, LR_para, LR_method_name, LR_method_filename)

        VFI_intercept = LR_para(:, 1);
        VFI_slope = LR_para(:, 2);
        RNFL_intercept = LR_para(:, 3);
        RNFL_slope = LR_para(:, 4);
                
        %% VFI intercept and RNFL intercept
        scatter(VFI_intercept,RNFL_intercept);
        xlabel('VFI intercept');
        ylabel('RNFL intercept');
        str = sprintf('%s (VFI intercept: %.1f +- %.1f, RNFL intercept: %.1f +- %.1f', LR_method_name, mean(VFI_intercept), std(VFI_intercept), mean(RNFL_intercept), std(RNFL_intercept));
        title(str);
        filename = sprintf('%s/%s_VFI_i_RNFL_i', out_dir, LR_method_filename);
        print(filename, '-dpng', '-r600'); %<-Save as PNG with 600 DPIs
        
        %% VFI slope and RNFL slope
        scatter(VFI_slope,RNFL_slope);
        xlabel('VFI slope');
        ylabel('RNFL slope');
        str = sprintf('%s (VFI slope: %.1f +- %.1f, RNFL slope: %.1f +- %.1f', LR_method_name, mean(VFI_slope), std(VFI_slope), mean(RNFL_slope), std(RNFL_slope));
        title(str);
        filename = sprintf('%s/%s_VFI_s_RNFL_s', out_dir, LR_method_filename);
        print(filename, '-dpng', '-r600'); %<-Save as PNG with 600 DPIs
        
        %% VFI intercept and slope
        scatter(VFI_slope, VFI_intercept);
        xlabel('VFI slope');
        ylabel('VFI intercept');
        str = sprintf('%s (VFI slope: %.1f +- %.1f, VFI intercept: %.1f +- %.1f', LR_method_name, mean(VFI_slope), std(VFI_slope), mean(VFI_intercept), std(VFI_intercept));
        title(str);
        filename = sprintf('%s/%s_VFI_s_VFI_i', out_dir, LR_method_filename);
        print(filename, '-dpng', '-r600'); %<-Save as PNG with 600 DPIs
        
        %% RNFL intercept and slope
        scatter(RNFL_slope, RNFL_intercept);
        xlabel('RNFL slope');
        ylabel('RNFL intercept');
        str = sprintf('%s (RNFL slope: %.1f +- %.1f, RNFL intercept: %.1f +- %.1f', LR_method_name, mean(RNFL_slope), std(RNFL_slope), mean(RNFL_intercept), std(RNFL_intercept));
        title(str);
        filename = sprintf('%s/%s_RNFL_s_RNFL_i', out_dir, LR_method_filename);
        print(filename, '-dpng', '-r600'); %<-Save as PNG with 600 DPIs
        
    end

end

