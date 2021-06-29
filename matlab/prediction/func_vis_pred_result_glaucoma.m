function func_vis_pred_result_glaucoma(subject_data, hist_state_seq, out_filename)

global data_setting;
global state_list;
global num_min_hist_visit;

type_name_ls = data_setting.type_name_ls;

figure,

set(gca, 'FontSize', 6);

if (num_min_hist_visit >= 2)
    global_LR_regress = subject_data.global_LR_regress;
    hist_LR_regress = subject_data.hist_LR_regress;
    hist_bayes_regress = subject_data.hist_bayes_regress;
end

ori_obs_time_seq = subject_data.ori_obs_time_seq;
ori_obs_seq = subject_data.ori_obs_seq;

ori_obs_time_seq = ori_obs_time_seq - ori_obs_time_seq(1);

num_hist_visit = size(hist_state_seq, 1);
num_dim = data_setting.dim;
CTHMM_pred_obs_seq = subject_data.CTHMM_pred_obs_seq;

for d = 1:num_dim
    
    %% set the subplot location
    subplot(num_dim, 1, d);
    
    set(gca, 'FontSize', 6);
    
    %a = get(gca,'XTickLabel');
    %set(gca,'XTickLabel',a,'fontsize',5);
    
    %a = get(gca,'YTickLabel');
    %set(gca,'YTickLabel',a,'fontsize',5);
    
            
    %% draw all original data points
    scatter(ori_obs_time_seq, ori_obs_seq(d, :), 90, 'filled', 'MarkerFaceColor', 'g'); % green solid dots  
    hold on;
      
    %% draw all CTHMM decoded states
    hist_HMM_decoded_obs = zeros(1, num_hist_visit);
    hist_HMM_half_range = zeros(1, num_hist_visit);
    
    for v = 1:num_hist_visit
        s = hist_state_seq(v);
        hist_HMM_decoded_obs(v) = state_list{s}.mu(d);
                
        dim_range_list = CTHMM_MD_query_dim_range_from_dim_idx(state_list{s}.dim_states);        
        hist_HMM_half_range(v) = abs(dim_range_list(d, 1) - dim_range_list(d, 2)) * 0.5;
    end
    
    X = ori_obs_time_seq(1:num_hist_visit);
    Y = hist_HMM_decoded_obs;
    E = hist_HMM_half_range;
    errorbar(X,Y,E, 'kd', 'LineWidth', 2); 
    
    %scatter(ori_obs_time_seq(1:num_hist_visit), hist_HMM_decoded_obs, 90, 'd', 'MarkerEdgeColor', 'k', 'LineWidth', 2); % black solid diamond
    hold on;
        
    %% draw CTHMM predicted data points
    predict_marker_color = [0.6 0.6 0.6];
    scatter(ori_obs_time_seq((num_hist_visit+1):end), CTHMM_pred_obs_seq(d, :)', 90, 'd', 'MarkerEdgeColor', predict_marker_color, 'LineWidth', 2); % gray diamond
    hold on; 
    
    
    if (num_min_hist_visit >= 2)
        
        %% draw regress line for the history
        hist_t_ls = [ori_obs_time_seq(1) ori_obs_time_seq(num_hist_visit)];
        future_t_ls = [ori_obs_time_seq(num_hist_visit) ori_obs_time_seq(end)];

        %% draw LR
        plot(hist_t_ls, hist_LR_regress(d,1) + hist_LR_regress(d, 2) * hist_t_ls, 'm', 'LineWidth', 2);
        plot(future_t_ls, hist_LR_regress(d,1) + hist_LR_regress(d, 2) * future_t_ls, '--m', 'LineWidth', 2);

        %% draw Bayes LR
        plot(hist_t_ls, hist_bayes_regress(d,1) + hist_bayes_regress(d, 2) * hist_t_ls, 'b', 'LineWidth', 2);
        plot(future_t_ls, hist_bayes_regress(d,1) + hist_bayes_regress(d, 2) * future_t_ls, '--b', 'LineWidth', 2);
    
    end
    
    %% draw global LR
    if (num_min_hist_visit >=2)
        plot(hist_t_ls, global_LR_regress(d,1) + global_LR_regress(d, 2) * hist_t_ls, '--c', 'LineWidth', 2);
        plot(future_t_ls, global_LR_regress(d,1) + global_LR_regress(d, 2) * future_t_ls, '--c', 'LineWidth', 2);
    end
    
    %% draw axis label and title        
    xlabel('Months', 'FontSize', 8);
    ylabel(type_name_ls{d},'FontSize', 8);
    
    if (num_min_hist_visit >=2)
        str = sprintf('ID: %d, %s, %s, age=%.1f [LR: slope:%.2f, err: %.1f][Bayes LR: slope:%.2f, err: %.1f][HMM: err:%.1f][global LR: slope:%.2f, err:%.1f]\n', ...            
                subject_data.ID,...
                subject_data.eye{1},...
                subject_data.visit_list{1}.Dx{1}, ...
                subject_data.visit_list{1}.age, ...
                hist_LR_regress(d, 2), subject_data.ave_LR_pred_err(d), ...
                hist_bayes_regress(d, 2), subject_data.ave_bayes_pred_err(d), ...
                subject_data.ave_CTHMM_pred_err(d), ...
                global_LR_regress(d, 2), subject_data.ave_global_pred_err(d));
    else
        str = sprintf('ID: %d, %s, %s, age=%.1f [HMM: err:%.1f]\n', ...            
            subject_data.ID,...
            subject_data.eye{1},...
            subject_data.visit_list{1}.Dx{1}, ...
            subject_data.visit_list{1}.age, ...            
            subject_data.ave_CTHMM_pred_err(d));            
    end
    
    title(str, 'FontSize', 5);        
    legend({'Original data', 'HMM decoded state (history)', 'HMM prediction (future)', 'Linear regression (history)', 'Linear regression (future)', 'Bayesian linear regression (history)', 'Bayesian linear regression (future)', 'Global linear regression (all data)'}, 'Location', 'northeastoutside', 'FontSize', 5);
        
end

%saveas(gcf, out_filename);
print(out_filename, '-dpng', '-r600'); %<-Save as PNG with 600 DPIs
close(gcf);
