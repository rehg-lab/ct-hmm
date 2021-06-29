function CTHMM_vis_2D_most_probable_future_path(top_out_folder, state_count_threshold, future_time)

global state_list;
global data_setting;
global Nij_mat;


%% quiver function usage
%% http://stackoverflow.com/questions/25729784/how-to-draw-an-arrow-in-matlab
%drawArrow = @(x,y,varargin) quiver( x(1),y(1),x(2)-x(1),y(2)-y(1),0, varargin{:} );

figure,
num_state = size(state_list, 1);

%% draw setting
axis_fontsize = 4;
x_pos_scale = 2.5;
y_pos_scale = 1.5;
x_pos_offset = 1;
y_pos_offset = -1.5;
max_marker_size = 10;
min_marker_size = 2;

%=========================================================================
%% draw all states (size: ni, color: dwell time)
Ni_list = diag(Nij_mat);
min_Ni = min(Ni_list);
max_Ni = max(Ni_list);

%=========================================================================
%% draw future state arrow from each state
%% compute Pt
%========================================================================
%% draw all states at once
vis_threshold = 0.000001;

%=========================================================================
%% draw future path
for m = 1:num_state
        
    Ni = Ni_list(m);
    
    %% if state m is a end-state, use the count from the incoming links   
    Ni_incoming = sum(Nij_mat(:, m)) - Ni;        
    if (Ni_incoming < state_count_threshold && Ni < state_count_threshold)
        continue;
    end    
        
    figure,
    
    m
        
    dim_range_list = CTHMM_MD_query_dim_range_from_dim_idx(state_list{m}.dim_states);
    
    %%%%%%%%%%%%%%%%%%%%%%
    %% draw all states 
    
    for k = 1:num_state
    
        Ni = Ni_list(k);    
        %% if state m is a end-state, use the count from the incoming links   
        Ni_incoming = sum(Nij_mat(:, k));        
        if (Ni_incoming < vis_threshold && Ni < vis_threshold)
            continue;
        end    
        % proportional based on data count        
        marker_size = min_marker_size + double(Ni - min_Ni) / double(max_Ni - min_Ni) * double(max_marker_size - min_marker_size);

        %% current state position
        k_state_x = state_list{k}.dim_states(1) * x_pos_scale + x_pos_offset;
        if (data_setting.draw_origin_lefttop == 1)        
            k_state_y = -state_list{k}.dim_states(2) * y_pos_scale + y_pos_offset;
        else        
            k_state_y = state_list{k}.dim_states(2) * y_pos_scale + y_pos_offset;
        end        
        
        if (k ~= m)
            plot(k_state_x,k_state_y,'yo','MarkerEdgeColor','k','MarkerFaceColor',[1.0 1.0 1.0],'MarkerSize',marker_size);        
        else
            plot(k_state_x,k_state_y,'yo','MarkerEdgeColor','k','MarkerFaceColor',[0.8 0.8 0.8],'MarkerSize',marker_size);        
        end
        
        hold on;
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%
    %% compute the best future path up to time T
    start_s = m;
    T = future_time;
    tic;
    str = 'Run CTMC_decode_most_probable_future_path_SSA()...'
    
    [best_state_seq_SSA, best_prob_SSA] = CTMC_decode_most_probable_future_path_SSA(start_s, future_time);    
    tEnd = toc;
    str = sprintf('SSA: %d min, %f sec\n', floor(tEnd/60),rem(tEnd,60))
    
    %% expm method
    %[dur_list] = CTHMM_decode_expected_dur_for_a_path_Expm(best_state_seq_SSA, T);    
    %%%%%%%%%%%%%%%%%%%%%%%
    
    %% now draw the predicted future path
    num_state_in_path = length(best_state_seq_SSA);
    
    %cur_time = 0.0;
    
    for i = 1:1:(num_state_in_path-1)
        
        %% pre state position
        cur_s = best_state_seq_SSA(i);
        next_s = best_state_seq_SSA(i+1);
        
        cur_state_x = state_list{cur_s}.dim_states(1) * x_pos_scale + x_pos_offset;
        if (data_setting.draw_origin_lefttop == 1)        
            cur_state_y = -state_list{cur_s}.dim_states(2) * y_pos_scale + y_pos_offset;
        else        
            cur_state_y = state_list{cur_s}.dim_states(2) * y_pos_scale + y_pos_offset;
        end
        next_state_x = state_list{next_s}.dim_states(1) * x_pos_scale + x_pos_offset;
        if (data_setting.draw_origin_lefttop == 1)        
            next_state_y = -state_list{next_s}.dim_states(2) * y_pos_scale + y_pos_offset;
        else        
            next_state_y = state_list{next_s}.dim_states(2) * y_pos_scale + y_pos_offset;
        end
        start_point = [cur_state_x cur_state_y];
        width = (next_state_x - cur_state_x);            
        height = (next_state_y - cur_state_y);            
        h=annotation('arrow');
        set(h,'parent', gca, ...
            'position', [start_point(1) start_point(2) width height], ...
            'HeadLength', 8, 'HeadWidth', 8, ...
            'Color', [0.4 0.1 0.8], 'LineWidth', 1);

        hold on;
        
    end
    
%     %% now draw dwell time for all states in path
%     for i = 1:1:num_state_in_path         
%         %% pre state position
%         cur_s = best_state_seq_SSA(i);        
%         cur_state_x = state_list{cur_s}.dim_states(1) * x_pos_scale + x_pos_offset;
%         if (data_setting.draw_origin_lefttop == 1)        
%             cur_state_y = -state_list{cur_s}.dim_states(2) * y_pos_scale + y_pos_offset;
%         else        
%             cur_state_y = state_list{cur_s}.dim_states(2) * y_pos_scale + y_pos_offset;
%         end
%     
%         %% draw dwelling time      
%         dwell_time = dur_list(i) / 12.0;
%         %cur_time = cur_time + dwell_time;
%         
%         %% write text: dwell time
%         if (dwell_time > 100)
%             str = sprintf('>100');
%         else
%             str = sprintf('%.2f', dwell_time);
%         end        
%         % write dwell time
%         fontsize = 5;
%         text(cur_state_x+0.2,cur_state_y-0.1, str, 'FontSize', fontsize, 'Color', [0 0 0]);       
%         hold on;
%     end
              
    %=========================================================================
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
        if ((data_setting.dim_value_range_ls{1}(i) - floor(data_setting.dim_value_range_ls{1}(i))) > 0)
            label = sprintf('[%.1f-%.1f]', data_setting.dim_value_range_ls{1}(i), data_setting.dim_value_range_ls{1}(i+1));
            if (i == 1)
                text(i * x_pos_scale + x_pos_offset - 0.35, -0.8 ,label, 'Color', color, 'FontSize',axis_fontsize);
            else
                text(i * x_pos_scale + x_pos_offset, -0.8 ,label, 'Color', color, 'FontSize',axis_fontsize);
            end
        else
            label = sprintf('[%d-%d]', data_setting.dim_value_range_ls{1}(i), data_setting.dim_value_range_ls{1}(i+1));
            text(i * x_pos_scale + x_pos_offset, -0.8, label, 'Color', color, 'FontSize',axis_fontsize);
        end
    end
    for i = 1:data_setting.dim_state_num_ls(2) % dim2:y

        label = sprintf('[%d-%d]', data_setting.dim_value_range_ls{2}(i), data_setting.dim_value_range_ls{2}(i+1));
        %label = sprintf('[%.2f-%.2f]', data_setting.dim_value_range_ls{2}(i), data_setting.dim_value_range_ls{2}(i+1));    
        if (data_setting.draw_origin_lefttop == 1)
            text(0.65, -i * y_pos_scale + y_pos_offset, label, 'Color', color, 'FontSize', axis_fontsize);    
        else
            text(0.65, i * y_pos_scale + y_pos_offset, label, 'Color', color, 'FontSize', axis_fontsize);   
        end
    end

    future_year = round(future_time / 12.0);
    str = sprintf('The Most Probable Future State Path In %d Years', future_year);
    title(str);

    %% save files
    filename = sprintf('%s/vis_futurepath_%d_y_%d_%d_%d_%d', top_out_folder, future_year, round(dim_range_list(1,1)), round(dim_range_list(1,2)), round(dim_range_list(2,1)), round(dim_range_list(2,2)));    
    print(filename, '-dpng', '-r600'); %<-Save as PNG with 600 DPIs

    close(gcf);

end

