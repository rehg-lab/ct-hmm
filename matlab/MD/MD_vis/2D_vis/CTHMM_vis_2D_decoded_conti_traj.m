function  CTHMM_vis_2D_decoded_conti_traj(subject_data, decoded_conti_state_seq, decoded_conti_dwelltime_seq, decoded_obs_state_seq, out_filename)

scrsz = get(groot,'ScreenSize');
figure('Position', [scrsz(3)/6 scrsz(4)/6  scrsz(3)*2.7/5 scrsz(4)*3/4]);

%% plot the decoded path
global state_list;
global data_setting;

%% draw setting
fontsize = 8;
axis_fontsize = 14;
max_line_thick = 5;

x_pos_scale = 1;
y_pos_scale = 1;
x_pos_offset = 0;
y_pos_offset = 0;

num_conti_state = size(decoded_conti_state_seq, 1);
marker_size = 12;

%% set the subplot location
position_vec = [0.11, 0.05, 0.50, 0.8];
subplot('Position', position_vec);
set(gca, 'FontSize', fontsize);

%% draw continuous state sequence
for i = 1:num_conti_state
                   
    m = decoded_conti_state_seq(i);
    
    %%==============================================================
    %% draw line between decoded states         
    if (i ~= num_conti_state) % not the last state
        
        n = decoded_conti_state_seq(i+1);
        
        m_states = state_list{m}.dim_states; % dimension index
        n_states = state_list{n}.dim_states;

        line_color = [0 0 0]; % draw black link
        line_pattern = 'k';
        line_width = max_line_thick-2;

        if (data_setting.draw_origin_lefttop == 1)
            plot([m_states(1) n_states(1)].* x_pos_scale + x_pos_offset, [-m_states(2) -n_states(2)].* y_pos_scale + y_pos_offset, line_pattern,  'LineWidth', line_width, 'color', line_color);
        else
            plot([m_states(1) n_states(1)].* x_pos_scale + x_pos_offset, [m_states(2) n_states(2)].* y_pos_scale + y_pos_offset, line_pattern,  'LineWidth', line_width, 'color', line_color);
        end            
        hold on;
    end    
    %================================================================
        
    %% draw decoded state with a large black circle
    if (data_setting.draw_origin_lefttop == 1)
        draw_i = state_list{m}.dim_states(1) * x_pos_scale + x_pos_offset;
        draw_j = -state_list{m}.dim_states(2) * y_pos_scale + y_pos_offset;
    else
        draw_i = state_list{m}.dim_states(1) * x_pos_scale + x_pos_offset;
        draw_j = state_list{m}.dim_states(2) * y_pos_scale + y_pos_offset;
    end        
    
    dwell_time = decoded_conti_dwelltime_seq(i) / 12.0; % month to year unit 
    color = CTHMM_vis_compute_dwell_time_color(dwell_time);
    
    plot(draw_i, draw_j, 'ro','MarkerEdgeColor', 'k', 'MarkerFaceColor', color, 'MarkerSize', marker_size);  %'none'
    hold on;
    
    %% write text: expected dwelling time    
    str = sprintf('%.1fy', dwell_time);
    text(draw_i+0.25,draw_j+0.15, str, 'FontSize', 9, 'Color', [0 0 0], 'Rotation', 50);
                       
end

%% draw discrete decoded states at observation times
num_obs = size(decoded_obs_state_seq, 1);
i = 1;
while (i <= num_obs)
    
    cur_state_idx = decoded_obs_state_seq(i);
    j = i;
    end_same_state_pos = i;
    while (1)
        if (decoded_obs_state_seq(j) == cur_state_idx)
            end_same_state_pos = j;
        else
            break;
        end
        j = j + 1;
        if (j > num_obs)
            break;
        end
    end
    
    if (i ~= end_same_state_pos)
        str = sprintf('v%d-%d', i, end_same_state_pos);
    else
        str = sprintf('v%d', i);
    end
    %% write visit info
    m = cur_state_idx;
    if (data_setting.draw_origin_lefttop == 1)
        draw_i = state_list{m}.dim_states(1) * x_pos_scale + x_pos_offset;
        draw_j = -state_list{m}.dim_states(2) * y_pos_scale + y_pos_offset;
    else
        draw_i = state_list{m}.dim_states(1) * x_pos_scale + x_pos_offset;
        draw_j = state_list{m}.dim_states(2) * y_pos_scale + y_pos_offset;
    end 
    
    text(draw_i+0.15,draw_j+0.53, str, 'FontSize', 9, 'Color', [0 0 1], 'Rotation', 50);
    
    %% update the next check pos
    i = end_same_state_pos + 1;
    
end

%% draw axis label
type_name_ls = data_setting.type_name_ls;
str = sprintf('Functional degeneration (%s data)', type_name_ls{1});
xlabel(str, 'FontSize', axis_fontsize);
str = sprintf('Structural degeneration (%s data)', type_name_ls{2});
ylabel(str, 'FontSize', axis_fontsize);
axis equal;

%title(train_group_name);

if (data_setting.draw_origin_lefttop == 1)
    xlim([0 (data_setting.dim_state_num_ls(1) * x_pos_scale + 1 + x_pos_offset * 2)]);
    ylim([(-data_setting.dim_state_num_ls(2) * y_pos_scale - 1 + y_pos_offset * 2) 0]);
else
    xlim([0 (data_setting.dim_state_num_ls(1) * x_pos_scale + 1 + x_pos_offset * 2)]);
    ylim([0 (data_setting.dim_state_num_ls(2) * y_pos_scale + 1 + y_pos_offset * 2)]);    
end

%% draw state definition
%color = [0.0 0.0 0.0];

ax = gca;
ax.XTickLabel = cell(data_setting.dim_state_num_ls(1), 1);
ax.YTickLabel = cell(data_setting.dim_state_num_ls(2), 1);

ax.XTick = [1:1:data_setting.dim_state_num_ls(1)];
%ax.YTick = [-1:-1:-data_setting.dim_state_num_ls(2)];
ax.YTick = [-data_setting.dim_state_num_ls(2):1:-1];

ax.XTickLabelRotation = 30;
ax.YTickLabelRotation = 30;

for i = 1:data_setting.dim_state_num_ls(1) % dim1:x

    if ((data_setting.dim_value_range_ls{1}(i) - floor(data_setting.dim_value_range_ls{1}(i))) > 0)
        label = sprintf('[%.1f-%.1f]', data_setting.dim_value_range_ls{1}(i), data_setting.dim_value_range_ls{1}(i+1));        
    else
        label = sprintf('[%d-%d]', data_setting.dim_value_range_ls{1}(i), data_setting.dim_value_range_ls{1}(i+1));        
    end
    ax.XTickLabel{i} = label;
end

y_idx = 1;
for i = data_setting.dim_state_num_ls(2):(-1):1 % dim2:y
    label = sprintf('[%d-%d]', data_setting.dim_value_range_ls{2}(i), data_setting.dim_value_range_ls{2}(i+1));    
    ax.YTickLabel{y_idx} = label;
    y_idx = y_idx + 1;
end

str = sprintf('Decoded State Trajectory (baseline age: %.1f)', subject_data.visit_list{1}.age);
title(str, 'FontSize', axis_fontsize);
set(gca,'xaxisLocation','top');

%% draw function and structural degeneration arrow

start_point = [0.22 -0.2];
width = data_setting.dim_state_num_ls(1);
height = 0.0;
h=annotation('arrow');
set(h,'parent', gca, ...
    'position', [start_point(1) start_point(2) width height], ...
    'HeadLength', 8, 'HeadWidth', 8, ...
    'Color', [0.0 0.0 0.0], 'LineWidth', 2);

%text(14.5, -0.3, 'Functional degeneration', 'FontSize', 9, 'Color', [0 0 1]);

start_point = [0.22 -0.2];
width = 0;
height = -data_setting.dim_state_num_ls(2);
h=annotation('arrow');
set(h,'parent', gca, ...
    'position', [start_point(1) start_point(2) width height], ...
    'HeadLength', 8, 'HeadWidth', 8, ...
    'Color', [0.0 0.0 0.0], 'LineWidth', 2);

hold on;

%========================================================================

%% set the subplot location
type_name_ls = data_setting.type_name_ls;

%% compute subject global linear regression parameters
begin_visit_idx = 1;
num_visit = subject_data.num_visit;
end_visit_idx = num_visit;
time_origin_idx = 1;

temp_subject_data = subject_data;
for v = 1:num_visit
    temp_subject_data.visit_list{v}.time = temp_subject_data.visit_list{v}.time / 12.0; % change from month to year unit, as we like to draw linear regression in year unit
end
[global_LR_regress, global_LR_err_sigma_vec] = func_compute_linear_regress_para(temp_subject_data, begin_visit_idx, end_visit_idx, time_origin_idx);    
    
ori_obs_time_seq = subject_data.visit_time_list / 12.0; % year unit
ori_obs_seq = subject_data.visit_data_list;
ori_obs_time_seq = ori_obs_time_seq - ori_obs_time_seq(1);

num_dim = data_setting.dim;

for d = 1:num_dim
    
    %% set the subplot location    
    if (d == 1)        
        position_vec = [0.70, 0.55, 0.25, 0.30];
        subplot('Position', position_vec);
    elseif (d == 2)        
        position_vec = [0.70, 0.10, 0.25, 0.30];
        subplot('Position', position_vec);
    end        
    hold on;
    
    %% draw all original data points
    scatter(ori_obs_time_seq, ori_obs_seq(:, d), 90, 'filled', 'MarkerFaceColor', 'b'); % green solid dots  
    hold on;
      
    %% draw all CTHMM decoded states
    HMM_decoded_obs = zeros(1, num_visit);
    HMM_half_range = zeros(1, num_visit);
    
    for v = 1:num_visit
        s = decoded_obs_state_seq(v);
        HMM_decoded_obs(v) = state_list{s}.mu(d);
                
        dim_range_list = CTHMM_MD_query_dim_range_from_dim_idx(state_list{s}.dim_states);        
        HMM_half_range(v) = abs(dim_range_list(d, 1) - dim_range_list(d, 2)) * 0.5;
    end
    
    X = ori_obs_time_seq(1:num_visit);
    Y = HMM_decoded_obs;
    E = HMM_half_range;
    errorbar(X,Y,E, 'kd', 'LineWidth', 2);
        
    hold on;
       
    %% draw global LR    
    %% compute regression on all visits    
    plot(ori_obs_time_seq, global_LR_regress(d,1) + global_LR_regress(d, 2) * ori_obs_time_seq, '--b', 'LineWidth', 2);
    
    %% draw axis label and title        
    xlabel('Year', 'FontSize', axis_fontsize);
    ylabel(type_name_ls{d}, 'FontSize', axis_fontsize);  
    
    %subject_data.visit_list{1}.Dx{1},
    str = sprintf('%s data (LR slope:%.2f)', type_name_ls{d}, global_LR_regress(d, 2));
    title(str, 'FontSize', axis_fontsize);    
    
    %legend({'raw data', 'decoded state', 'LR slope'}, 'Location', 'northeastoutside');
    legend({'raw data', 'decoded state', 'LR line'}, 'Location', 'Best');
    
end

%========================================================================
%% save files
saveas(gcf, out_filename, 'png');
close(gcf);
