function func_vis_2D_raw_data_hist(top_out_dir, dim1_edges, dim2_edges, raw_data_2D_hist)


data_setting.draw_origin_lefttop = 1;

figure,

%% draw setting
fontsize = 7;
axis_fontsize = 4;

max_marker_size = 10;
min_marker_size = 1;

x_pos_scale = 2.5;
y_pos_scale = 1.5;
x_pos_offset = 1;
y_pos_offset = -1.5;

min_Si = min(raw_data_2D_hist(:));
max_Si = max(raw_data_2D_hist(:));

num_dim1_bin = length(dim1_edges) - 1;
num_dim2_bin = length(dim2_edges) - 1;

for x = 1:num_dim1_bin
    
    for y = 1:num_dim2_bin

        Si = raw_data_2D_hist(x, y);
        
        if (Si == 0) 
            continue;
        end
        
        % proportional based on data count        
        marker_size = min_marker_size + double(Si - min_Si) / double(max_Si - min_Si) * double(max_marker_size - min_marker_size);
    
        % draw state position
        if (data_setting.draw_origin_lefttop == 1)
            draw_i = x * x_pos_scale + x_pos_offset;
            draw_j = -y * y_pos_scale + y_pos_offset;
        else
            draw_i = x * x_pos_scale + x_pos_offset;
            draw_j = -y * y_pos_scale + y_pos_offset;
        end
    
        plot(draw_i, draw_j, 'ro','MarkerEdgeColor', 'k', 'MarkerFaceColor', 'y', 'MarkerSize', marker_size);
        label = sprintf('%d', Si); % raw data count
        str_color = [0 0 1];
        text(draw_i+0.2, draw_j+0.25, label, 'Color', str_color, 'FontSize', fontsize);        
        hold on;
        
    end    
end


%% draw axis label
%type_name_ls = data_setting.type_name_ls;
%xlabel(type_name_ls{1});
%ylabel(type_name_ls{2});
axis equal;


if (data_setting.draw_origin_lefttop == 1)
    xlim([0 (num_dim1_bin * x_pos_scale + 1 + x_pos_offset * 2)]);
    ylim([(-num_dim2_bin * y_pos_scale - 1 + y_pos_offset * 2) 0]);
else
    xlim([0 (num_dim1_bin * x_pos_scale + 1 + x_pos_offset * 2)]);
    ylim([0 (num_dim2_bin * y_pos_scale + 1 + y_pos_offset * 2)]);    
end

%% draw state definition
color = [0.0 0.0 0.0];
for i = 1:num_dim1_bin % dim1:x
    if ((dim1_edges(i) - floor(dim1_edges(i))) > 0)
        label = sprintf('[%.1f-%.1f]', dim1_edges(i), dim1_edges(i+1));
        if (i == 1)
            text(i * x_pos_scale + x_pos_offset - 0.35, -0.8 ,label, 'Color', color, 'FontSize',axis_fontsize);
        else
            text(i * x_pos_scale + x_pos_offset, -0.8 ,label, 'Color', color, 'FontSize',axis_fontsize);
        end
    else
        label = sprintf('[%d-%d]', dim1_edges(i), dim1_edges(i+1));
        text(i * x_pos_scale + x_pos_offset, -0.8, label, 'Color', color, 'FontSize',axis_fontsize);
    end
end
for i = 1:num_dim2_bin % dim2:y    
    label = sprintf('[%d-%d]',dim2_edges(i), dim2_edges(i+1));
    %label = sprintf('[%.2f-%.2f]', data_setting.dim_value_range_ls{2}(i), data_setting.dim_value_range_ls{2}(i+1));    
    if (data_setting.draw_origin_lefttop == 1)
        text(0.65, -i * y_pos_scale + y_pos_offset, label, 'Color', color, 'FontSize', axis_fontsize);    
    else
        text(0.65, i * y_pos_scale + y_pos_offset, label, 'Color', color, 'FontSize', axis_fontsize);   
    end
end

%% save files
filename = sprintf('%s/vis_2D_raw_data_hist', top_out_dir);
print(filename, '-dpng', '-r600'); %<-Save as PNG with 600 DPIs
%close(gcf);

