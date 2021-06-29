function CTHMM_vis_2D_most_probable_future_state(top_out_folder, vis_threshold, future_time)

global state_list;
global Q_mat;
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
Pt = expm(future_time * Q_mat);

%========================================================================
%% draw all states dot at once
for m = 1:num_state
    
    Ni = Ni_list(m);    
    %% if state m is a end-state, use the count from the incoming links   
    Ni_incoming = sum(Nij_mat(:, m));        
    if (Ni_incoming < vis_threshold && Ni < vis_threshold)
        continue;
    end    
    % proportional based on data count        
    marker_size = min_marker_size + double(Ni - min_Ni) / double(max_Ni - min_Ni) * double(max_marker_size - min_marker_size);
    
    %% current state position
    m_state_x = state_list{m}.dim_states(1) * x_pos_scale + x_pos_offset;
    if (data_setting.draw_origin_lefttop == 1)        
        m_state_y = -state_list{m}.dim_states(2) * y_pos_scale + y_pos_offset;
    else        
        m_state_y = state_list{m}.dim_states(2) * y_pos_scale + y_pos_offset;
    end
    
    %% draw current state     
    temp = Pt(m, :);            
    [C, n] = max(temp);    
    if (n == m)
        marker_color = 'g';
    else
        marker_color = 'm';
    end
    plot(m_state_x,m_state_y,'yo','MarkerEdgeColor','k','MarkerFaceColor',marker_color,'MarkerSize',marker_size);        
    hold on;       
end

%=========================================================================
%% draw all future state at once
for m = 1:num_state
        
    Ni = Ni_list(m);
    
    %% if state m is a end-state, use the count from the incoming links   
    Ni_incoming = sum(Nij_mat(:, m));        
    if (Ni_incoming < vis_threshold && Ni < vis_threshold)
        continue;
    end    
        
    %% current state position
    m_state_x = state_list{m}.dim_states(1) * x_pos_scale + x_pos_offset;
    if (data_setting.draw_origin_lefttop == 1)        
        m_state_y = -state_list{m}.dim_states(2) * y_pos_scale + y_pos_offset;
    else        
        m_state_y = state_list{m}.dim_states(2) * y_pos_scale + y_pos_offset;
    end
        
    temp = Pt(m, :);            
    [C, n] = max(temp);
    
    %% find future state n's position        
    if (n~= m && C ~= 0)
        I = find(temp == C);
        num_strong_link = size(I, 2);
        
        for k = 1:num_strong_link
            n = I(k);                
            n_state_x = state_list{n}.dim_states(1) * x_pos_scale + x_pos_offset;
            if (data_setting.draw_origin_lefttop == 1)                        
                n_state_y = -state_list{n}.dim_states(2) * y_pos_scale + y_pos_offset;
            else                    
                n_state_y = -state_list{n}.dim_states(2) * y_pos_scale + y_pos_offset;
            end                      
            start_point = [m_state_x m_state_y];
            width = (n_state_x - m_state_x);            
            height = (n_state_y-m_state_y);
            
            h=annotation('arrow');
            set(h,'parent', gca, ...
                'position', [start_point(1) start_point(2) width height], ...
                'HeadLength', 8, 'HeadWidth', 8, ...
                'Color', [0.4 0.1 0.8], 'LineWidth', 1);
                                    
            hold on;
        end                        
    end 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%        
end

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
str = sprintf('The Most Probable Future State After %d Years', future_year);
title(str);


%% save files
filename = sprintf('%s/vis_futurestate_%d_y', top_out_folder, future_year);
%print('-dtiff','-r600',filename);
print(filename, '-dpng', '-r600'); %<-Save as PNG with 600 DPIs

close(gcf);

