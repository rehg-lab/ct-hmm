function CTHMM_vis_2D_Q_mat_new(top_out_folder, vis_mat_type, vis_threshold, is_draw_text, state_is_vis_tau_i)

% state_is_vis_tau_i = 1, we draw state size based on tau_i
% state_is_vis_tau_i = 0, we draw state size based on n_i = \sum n_{ij}

global state_list;
global Q_mat_struct;
global Q_mat;
global data_setting;
%global train_group_name;

global Nij_mat;
global Ti_list;


num_nij_less_neg001 = length(find(Nij_mat(:) < -0.01));

num_nij_less_neg001

figure,

%% draw setting
%fontsize = 3;
%axis_fontsize = 4;

fontsize = 2;
axis_fontsize = 2;


%max_marker_size = 14;
%min_marker_size = 2;
%max_line_thick = 5;
%min_line_thick = 0.25;

max_marker_size = 10;
%min_marker_size = 2;
min_marker_size = 1;
max_line_thick = 3;
min_line_thick = 0.25;

%x_pos_scale = 2.5;
%y_pos_scale = 1.5;

x_pos_scale = 3;
y_pos_scale = 2;

%x_pos_offset = 1;
x_pos_offset = 6;
y_pos_offset = -2;

%% plot transition intensity (qrs) at each arrow
num_state = size(state_list, 1);

if (vis_mat_type == 1) % vis nij
    vis_mat = Nij_mat;
    max_line_thick = 3;
    min_line_thick = 0.25;
elseif (vis_mat_type == 2) % vis qij     
    vis_mat = Q_mat;
    max_line_thick = 3;
    min_line_thick = 0.25;
elseif (vis_mat_type == 3) % vis vij
    vis_mat = Nij_mat;
    for s = 1:num_state
        vis_mat(s,s) = 0;
        sum_row = sum(vis_mat(s, :));
        if (sum_row ~= 0)
            vis_mat(s, :) = vis_mat(s, :) ./ sum_row;
        end
    end
    max_line_thick = 1.5;
    min_line_thick = 0.25;
end

%% find maximum value in the vis mat
temp_mat = vis_mat;
for s = 1:num_state
    temp_mat(s,s) = 0;
end
max_vis_ij = max(temp_mat(:));

%% plot nij or qij link
for m = 1:num_state
    
    sum_link = sum(Q_mat_struct(m, :));    
    if (sum_link == 0)
        continue;
    end
        
    vis_row = vis_mat(m, :);
    vis_row(1, m) = 0;    
    [C, temp] = max(vis_row); % max value in vis row
        
    for n = 1:num_state
        
        if (Q_mat_struct(m, n) == 1)
            
            if (vis_mat_type == 1 && vis_mat(m, n) < vis_threshold)
                continue;
            end
            
            % dimension index
            m_states = state_list{m}.dim_states;
            n_states = state_list{n}.dim_states;
            
            %% draw link position
%            x = (m_states(1) + n_states(1))/2.0 * x_pos_scale + x_pos_offset;
%            y = (m_states(2) + n_states(2))/2.0 * y_pos_scale + y_pos_offset;                
            
            infov = vis_mat(m, n);

            %% draw links
            if (vis_mat(m, n) == C && C > 0)  % the strongest link
                line_color = [0 0 1]; % draw blue link
                line_pattern = 'b';
            else
                line_color = [0 0 0]; % draw black link
                line_pattern = '--k';
            end
            
            if (vis_mat_type == 1 || vis_mat_type == 2)
                line_width = min_line_thick + double(infov) / double(max_vis_ij) * double(max_line_thick - min_line_thick);
            else
                line_width = min_line_thick + infov * double(max_line_thick - min_line_thick);
            end
                                    
            if (data_setting.draw_origin_lefttop == 1)
                plot([m_states(1) n_states(1)].* x_pos_scale + x_pos_offset, [-m_states(2) -n_states(2)].* y_pos_scale + y_pos_offset, line_pattern,  'LineWidth', line_width, 'color', line_color);
            else
                plot([m_states(1) n_states(1)].* x_pos_scale + x_pos_offset, [m_states(2) n_states(2)].* y_pos_scale + y_pos_offset, line_pattern,  'LineWidth', line_width, 'color', line_color);
            end            
            hold on;                        
        end % if
    end % n
end

%% draw all states
%

if (state_is_vis_tau_i == 1)
    state_vis_list = Ti_list;
else
    Ni_list = diag(Nij_mat);
    state_vis_list = Ni_list;
end

min_Si = min(state_vis_list);
max_Si = max(state_vis_list);

for m = 1:num_state
        
    Si = state_vis_list(m);      
    
    % proportional based on data count        
    marker_size = min_marker_size + double(Si - min_Si) / double(max_Si - min_Si) * double(max_marker_size - min_marker_size);
    
    % draw state position
    if (data_setting.draw_origin_lefttop == 1)
        draw_i = state_list{m}.dim_states(1) * x_pos_scale + x_pos_offset;
        draw_j = -state_list{m}.dim_states(2) * y_pos_scale + y_pos_offset;
    else
        draw_i = state_list{m}.dim_states(1) * x_pos_scale + x_pos_offset;
        draw_j = state_list{m}.dim_states(2) * y_pos_scale + y_pos_offset;
    end
    
    % plot average dwelling time
    ave_dwell_time = 1 / -Q_mat(m, m); % month
    ave_dwell_time = ave_dwell_time / 12.0; % month to year unit
    
    %% plot ave dwelling time color and state size
    %if (ave_dwell_time == inf) % absorping state
        %plot(draw_i, draw_j, 'ro', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [1.0 1.0 1.0], 'MarkerSize', marker_size);   
    %    plot(draw_i, draw_j, 'ro', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [0.0 0.0 0.0], 'MarkerSize', marker_size);   
    %else
        color = CTHMM_vis_compute_dwell_time_color(ave_dwell_time);
        plot(draw_i, draw_j, 'ro','MarkerEdgeColor', 'k', 'MarkerFaceColor', color, 'MarkerSize', marker_size);
    %end
    
    if (is_draw_text == 1)    
        %% write text: ave dwell time
        if (ave_dwell_time > 100)
            str = sprintf('>100y');
        else
            str = sprintf('%.2fy', ave_dwell_time);
        end        
        % write ave dwell time
        text(draw_i+0.2,draw_j-0.1, str, 'FontSize', fontsize, 'Color', [1 0 0]);

        %% write text: total outflow subject count (n_i) or expected total dwelling time at this state
        if (state_is_vis_tau_i == 1)            
            label = sprintf('%.2fy', Si / 12.0); % month to year unit
            str_color = [0.7 0.7 0];
        else
            label = sprintf('#%.1f', Si); % total outflow subject count (n_i)
            str_color = [0 0.5 0];
        end
        text(draw_i+0.2, draw_j+0.25, label, 'Color', str_color, 'FontSize', fontsize);
        
        hold on;
    end
            
end

%% draw link texts
if (is_draw_text == 1)  
    
%% plot nij or qij link
for m = 1:num_state
    
    sum_link = sum(Q_mat_struct(m, :));    
    if (sum_link == 0)
        continue;
    end
        
    vis_row = vis_mat(m, :);
    vis_row(1, m) = 0;    
    [C, temp] = max(vis_row); % max value in vis row
        
    for n = 1:num_state
        
        if (Q_mat_struct(m, n) == 1)
            
            if (vis_mat_type == 1 && vis_mat(m, n) < vis_threshold)
                continue;
            end
            
            % dimension index
            m_states = state_list{m}.dim_states;
            n_states = state_list{n}.dim_states;
            
            %% draw link position
            x = (m_states(1) + n_states(1))/2.0 * x_pos_scale + x_pos_offset;
                        
            infov = vis_mat(m, n);
            
            if (vis_mat_type == 1) % expected count
                draw_info = sprintf('#%.1f', infov);
            else % transition rate or prob
                draw_info = sprintf('%.2f', infov);  
            end  
                        
            %% draw info text                                  
            if (data_setting.draw_origin_lefttop == 1)
                y = (-m_states(2) -n_states(2))/2.0 * y_pos_scale + y_pos_offset;                
            else
                y = (m_states(2) + n_states(2))/2.0 * y_pos_scale + y_pos_offset;                
            end
            %text(x+0.1, y+0.15, draw_info, 'Color', [0 0.5 0], 'FontSize', fontsize);
            text(x+0.1, y+0.25, draw_info, 'Color', [0 0.5 0], 'FontSize', fontsize);
            
            hold on;
                        
        end % if
    end % n
end
end

%% draw axis label
type_name_ls = data_setting.type_name_ls;
xlabel(type_name_ls{1});
ylabel(type_name_ls{2});
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
color = [0.0 0.0 0.0];
for i = 1:data_setting.dim_state_num_ls(1) % dim1:x
    if (((data_setting.dim_value_range_ls{1}(i) - floor(data_setting.dim_value_range_ls{1}(i))) > 0) || ...
        ((data_setting.dim_value_range_ls{1}(i+1) - floor(data_setting.dim_value_range_ls{1}(i+1))) > 0))
        label = sprintf('[%.1f-%.1f]', data_setting.dim_value_range_ls{1}(i), data_setting.dim_value_range_ls{1}(i+1)); % if range is in floating point
        if (i == 1)
            text(i * x_pos_scale + x_pos_offset - 0.35, -2 ,label, 'Color', color, 'FontSize',axis_fontsize);
        else
            text(i * x_pos_scale + x_pos_offset, -2 ,label, 'Color', color, 'FontSize',axis_fontsize);
        end
    else
        label = sprintf('[%d-%d]', data_setting.dim_value_range_ls{1}(i), data_setting.dim_value_range_ls{1}(i+1));
        text(i * x_pos_scale + x_pos_offset, -2, label, 'Color', color, 'FontSize',axis_fontsize);
    end
end
for i = 1:data_setting.dim_state_num_ls(2) % dim2:y    
    
    if (((data_setting.dim_value_range_ls{2}(i) - floor(data_setting.dim_value_range_ls{2}(i))) > 0) || ...   % if range is in floating point
        ((data_setting.dim_value_range_ls{2}(i+1) - floor(data_setting.dim_value_range_ls{2}(i+1))) > 0))    
        label = sprintf('[%.1f-%.1f]', data_setting.dim_value_range_ls{2}(i), data_setting.dim_value_range_ls{2}(i+1));
    else
        label = sprintf('[%d-%d]', data_setting.dim_value_range_ls{2}(i), data_setting.dim_value_range_ls{2}(i+1));
    end        
    if (data_setting.draw_origin_lefttop == 1)
        text(0.65, -i * y_pos_scale + y_pos_offset, label, 'Color', color, 'FontSize', axis_fontsize);    
    else
        text(0.65, i * y_pos_scale + y_pos_offset, label, 'Color', color, 'FontSize', axis_fontsize);   
    end
    
end

%% save files

if (vis_mat_type == 1) % vis nij
    link_size_str = 'nij';
elseif (vis_mat_type == 2) % vis qij     
    link_size_str = 'qij';
elseif (vis_mat_type == 3) % vis vij
    link_size_str = 'vij';
end
    
if (state_is_vis_tau_i == 1)
    state_size_str = 'taui';
else
    state_size_str = 'ni';
end

filename = sprintf('%s/vis_%s_%s', top_out_folder, link_size_str, state_size_str);

%saveas(gcf, filename, 'png');
%saveas(gcf, filename, 'epsc');
%saved_fig = sprintf('%s.fig', filename);
%savefig(saved_fig);
%print('-dtiff','-r600',filename);
print(filename, '-dpng', '-r600'); %<-Save as PNG with 600 DPIs
close(gcf);

