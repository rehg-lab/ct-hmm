function func_vis_glaucoma_2D_Q_mat_colormap(top_out_folder, future_time_list, num_total_subject, is_overlay_text, min_count_to_draw_state, min_accum_dwell_time_to_draw_absorb_state)

% colormap: for each state, compute an HSV color value as below, and then
% tranform to RGB value to draw
% Hue: represent the mean direction of transitions (0-90 degree for the glaucoma
% model, where 0: horizontal, 45: diagonal, 90: vertical)
% Saturation: set 1
% Value: the value level represents the probability of exiting a particular
% state after the specified time (months) for normalization setting 1;
% alternatively, it represents the number of subjects exiting a state
% normalizaed by the entire population after the specified time

global Q_mat_struct;
global Q_mat;
global Nij_mat;

global state_list;
global data_setting;
global Ti_list;

num_state = length(state_list);

%% draw setting
fontsize = 2;
axis_fontsize = 2;

%% draw position setting
x_pos_scale = 3;
y_pos_scale = x_pos_scale;
x_pos_offset = 4;
y_pos_offset = -1;

%% computer Vij matrix from Nij
% Vij_mat : the probability of transition to the neighbors
Vij_mat = Nij_mat;
for s = 1:num_state
    Vij_mat(s,s) = 0;
    sum_row = sum(Vij_mat(s, :));
    if (sum_row ~= 0)
        Vij_mat(s, :) = Vij_mat(s, :) ./ sum_row;
    end
end

Ni = diag(Nij_mat);
max_ni = max(Ni);

%% create two output folders
out_folder = sprintf('%s/norm_by_state', top_out_folder);
mkdir(out_folder);
out_folder = sprintf('%s/norm_by_popu', top_out_folder);
mkdir(out_folder);

num_future_time = length(future_time_list);

%% for each state, compute the probability of exiting to the other states after a specified time,
%% and compute the mean direction of transitions

for time_idx = 1:num_future_time
        
    future_time = future_time_list(time_idx);
    
    %% compute the transition probability matrix
    Pt = expm(Q_mat * future_time);
        
    for norm_setting = 1:2  % normalization by state (1) or by the entire population (2)
        
        figure,
        if (norm_setting == 1)
            title_str = sprintf('Transition Probability Map After %.1f Years, Norm by Individual State', future_time/12.0);
        else
            title_str = sprintf('Transition Count/Percent Map After %.1f Years, Norm by Population', future_time/12.0);
        end
        
        for m = 1:num_state % for each state
            
            m_states = state_list{m}.dim_states;
            
            %% compute the state drawing location                        
            draw_x = state_list{m}.dim_states(1) * x_pos_scale + x_pos_offset;                
            if (data_setting.draw_origin_lefttop == 1)                
                draw_y = -state_list{m}.dim_states(2) * y_pos_scale + y_pos_offset;                
            else                
                draw_y = state_list{m}.dim_states(2) * y_pos_scale + y_pos_offset;                
            end
            rect = [draw_x - 0.5*x_pos_scale, draw_y - 0.5*y_pos_scale, x_pos_scale, y_pos_scale];    
            
            sum_link = sum(Q_mat_struct(m, :));    
            if (sum_link == 0)  % no outgoing link, this is an absorption state, let's draw it as a black block
                
                if (Ti_list(m) >= min_accum_dwell_time_to_draw_absorb_state)
                    %% plot the HSV to RGB value for this state
                    rgb_color = [0 0 0]; % black
                    rectangle('Position',rect,'FaceColor',rgb_color,'EdgeColor','none');
                end
                
                continue;                
            end

            %% if a state's count is less than the threshold, not drawing it
            if (Nij_mat(m,m) < min_count_to_draw_state)
                continue;
            end
                        
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% now compute mean direction
            neighbor_index_ls = find(Q_mat_struct(m, :) == 1);    
            num_neighbor = length(neighbor_index_ls);
            mean_degree = 0;
            
            %Vij_mat(m, neighbor_index_ls)
                        
            % for each neighbor, accumulate the transition direction by
            % probability
            for i = 1:num_neighbor

                %% get neighbor's state index
                n = neighbor_index_ls(i);

                %% dimension index for the neighboring state        
                n_states = state_list{n}.dim_states;

                %% now check the neighbor state's direction, there are three possible cases for our 2D glaucoma progressing models                        
                if ((n_states(1) == (m_states(1) + 1)) && (n_states(2) == m_states(2)))  % n is the horizontal neighbor of m
                    %mean_degree = mean_degree + Vij_mat(m, n) * 0;
                elseif ((n_states(1) == (m_states(1) + 1)) && (n_states(2) == (m_states(2) + 1)))  % n is the diagonal neighbor of m
                    mean_degree = mean_degree + Vij_mat(m, n) * 45.0;
                elseif ((n_states(1) == m_states(1)) && (n_states(2) == (m_states(2) + 1))) % n is the vertical neighbor of m
                    mean_degree = mean_degree + Vij_mat(m, n) * 90.0;
                else 
                    % error 
                    disp('error in m,n state index');
                end                
            end

            %% find the transition probability for future time t from Pt
            prob_of_stay_m = Pt(m, m);
            prob_of_exiting_m = 1.0 - prob_of_stay_m;

            %% now assign the HSV value for state m
            if (norm_setting == 1) % normalization by state               
                value = prob_of_exiting_m;                
            else  % normalization by entire population 
                count_at_m = Nij_mat(m, m);
                %value = prob_of_exiting_m * count_at_m / double(num_total_subject);                
                value = prob_of_exiting_m * count_at_m / double(max_ni); % max ni is the maximum size at a state
            end
                                   
            HSV_color = [mean_degree/360.0 1.0 value];
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%           

            %% plot the HSV to RGB value for this state
            rgb_color = hsv2rgb(HSV_color);            
            rectangle('Position',rect,'FaceColor',rgb_color,'EdgeColor','none');
            
            if (is_overlay_text == 1)       
                                
                if (norm_setting == 1)
                    str = sprintf('%.1f%%\n%.1f%c', value*100, mean_degree, char(176));   
                else
                    count_out_of_m = prob_of_exiting_m * count_at_m;
                    percent_out_of_m = count_out_of_m / double(num_total_subject);
                    str = sprintf('#%.1f (%.1f%%)\n%.1f%c', count_out_of_m, percent_out_of_m*100, mean_degree, char(176));   
                end
                 
                if (value < 0.5)                
                    text(draw_x-0.30*x_pos_scale, draw_y, str, 'FontSize', fontsize, 'Color', [1 1 1]);
                else
                    text(draw_x-0.30*x_pos_scale, draw_y, str, 'FontSize', fontsize, 'Color', [0 0 0]);       
                end
                hold on;
            end

        end % for m

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
            else
                label = sprintf('[%d-%d]', data_setting.dim_value_range_ls{1}(i), data_setting.dim_value_range_ls{1}(i+1));
            end
            if (i == 1)
                text((i-0.30) * x_pos_scale + x_pos_offset - 0.4, -0.8 ,label, 'Color', color, 'FontSize',axis_fontsize);
            else
                text((i-0.30) * x_pos_scale + x_pos_offset, -0.8 ,label, 'Color', color, 'FontSize',axis_fontsize);
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
        
        %% save figures
        if (norm_setting == 1)
            out_folder = sprintf('%s/norm_by_state', top_out_folder);
            norm_str = 'normstate';
        elseif (norm_setting == 2)
            out_folder = sprintf('%s/norm_by_popu', top_out_folder);
            norm_str = 'normpopu';
        end
                       
        title(title_str);
        
        filename = sprintf('%s/vis_colormap_%dm_%s', out_folder, future_time, norm_str);    
        print(filename, '-dpng', '-r600'); %<-Save as PNG with 600 DPIs
        close(gcf);

     end % for norm_setting

end % for each future time
