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

        fontsize = 8;
        
        VFI_intercept = LR_para(:, 1);
        VFI_slope = LR_para(:, 2);
        RNFL_intercept = LR_para(:, 3);
        RNFL_slope = LR_para(:, 4);
                
        %% VFI intercept and RNFL intercept
        figure,
        scatter(VFI_intercept,RNFL_intercept);
        xlabel('VFI intercept');
        ylabel('RNFL intercept');
        str = sprintf('%s (VFI intercept: %.2f +- %.2f, RNFL intercept: %.2f +- %.2f)', LR_method_name, mean(VFI_intercept), std(VFI_intercept), mean(RNFL_intercept), std(RNFL_intercept));
        title(str, 'FontSize', fontsize);
        filename = sprintf('%s/%s_VFI_i_RNFL_i', out_dir, LR_method_filename);
        print(filename, '-dpng', '-r600'); %<-Save as PNG with 600 DPIs
        close(gcf);
        
        %% VFI slope and RNFL slope
        figure,
        scatter(VFI_slope,RNFL_slope);
        xlabel('VFI slope');
        ylabel('RNFL slope');
        str = sprintf('%s (VFI slope: %.2f +- %.2f, RNFL slope: %.2f +- %.2f)', LR_method_name, mean(VFI_slope), std(VFI_slope), mean(RNFL_slope), std(RNFL_slope));
        title(str, 'FontSize', fontsize);
        filename = sprintf('%s/%s_VFI_s_RNFL_s', out_dir, LR_method_filename);
        print(filename, '-dpng', '-r600'); %<-Save as PNG with 600 DPIs
        close(gcf);
        
        %% VFI intercept and slope
        figure,
        scatter(VFI_slope, VFI_intercept);
        xlabel('VFI slope');
        ylabel('VFI intercept');
        str = sprintf('%s (VFI slope: %.2f +- %.2f, VFI intercept: %.2f +- %.2f)', LR_method_name, mean(VFI_slope), std(VFI_slope), mean(VFI_intercept), std(VFI_intercept));
        title(str, 'FontSize', fontsize);
        filename = sprintf('%s/%s_VFI_s_VFI_i', out_dir, LR_method_filename);
        print(filename, '-dpng', '-r600'); %<-Save as PNG with 600 DPIs
        close(gcf);
        
        %% RNFL intercept and slope
        figure,
        scatter(RNFL_slope, RNFL_intercept);
        xlabel('RNFL slope');
        ylabel('RNFL intercept');
        str = sprintf('%s (RNFL slope: %.2f +- %.2f, RNFL intercept: %.2f +- %.2f)', LR_method_name, mean(RNFL_slope), std(RNFL_slope), mean(RNFL_intercept), std(RNFL_intercept));
        title(str, 'FontSize', fontsize);
        filename = sprintf('%s/%s_RNFL_s_RNFL_i', out_dir, LR_method_filename);
        print(filename, '-dpng', '-r600'); %<-Save as PNG with 600 DPIs
        close(gcf);
        
    end

end

