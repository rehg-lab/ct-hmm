function utility_plot_datasheet_statistics(fp_log, output_folder, output_prefix, num_visit_list, trace_year_list, VFI_list, RNFL_list, rimarea_list, baseline_age_list, age_list)

figure,
max_value = max(num_visit_list(:));
hist(num_visit_list, [1:1:max_value]);
str = sprintf('histogram of visit count (mean:%.1f, std=%.1f)', mean(num_visit_list), std(num_visit_list));
title(str);
str = sprintf('hist_visit_count_%s', output_prefix);
utility_save_high_quality_fig(gcf, output_folder, str);

figure,
max_value = max(trace_year_list(:));
hist(trace_year_list, [1:1:max_value]);
str = sprintf('histogram of tracing years(mean:%.1f, std=%.1f)', mean(trace_year_list), std(trace_year_list));
title(str);
str = sprintf('hist_tracing_year_%s', output_prefix);
utility_save_high_quality_fig(gcf, output_folder, str);

figure,
hist(VFI_list, [0:1:100]);
title('histogram of VFI');
xlim([0 100]);
str = sprintf('hist_VFI_%s', output_prefix);
utility_save_high_quality_fig(gcf, output_folder, str);

num_VFI_data = size(VFI_list, 2);
fprintf(fp_log, 'num VFI data = %d\n', num_VFI_data);

fprintf(fp_log, 'VFI == 100 (%.2f)\n', size(find(VFI_list == 100), 2)/double(num_VFI_data));
fprintf(fp_log, 'VFI ==  99 (%.2f)\n', size(find(VFI_list == 99), 2)/double(num_VFI_data));
fprintf(fp_log, 'VFI ==  98 (%.2f)\n', size(find(VFI_list == 98), 2)/double(num_VFI_data));
fprintf(fp_log, 'VFI >=  95 (%.2f)\n', size(find(VFI_list >= 95), 2)/double(num_VFI_data));
fprintf(fp_log, 'VFI >=  90 (%.2f)\n', size(find(VFI_list >= 90), 2)/double(num_VFI_data));

figure,
hist(RNFL_list, [0:5:155]);
title('histogram of RNFL');
xlim([0 155]);
fprintf(fp_log, 'min RNFL (%.1f)\n', min(RNFL_list));
fprintf(fp_log, 'max RNFL (%.1f)\n', max(RNFL_list));
str = sprintf('hist_RNFL_%s', output_prefix);
utility_save_high_quality_fig(gcf, output_folder, str);



figure,
hist(rimarea_list, [0.5:0.125:2.75]);
title('histogram of rim area');
xlim([0 3]);
fprintf(fp_log, 'min rimarea (%.1f)\n', min(rimarea_list));
fprintf(fp_log, 'max rimarea (%.1f)\n', max(rimarea_list));
str = sprintf('hist_rimarea_%s', output_prefix);
utility_save_high_quality_fig(gcf, output_folder, str);


figure,
hist(age_list, [0:10:100]);
str = sprintf('histogram of visit ages (mean:%.1f, std=%.1f)', mean(age_list), std(age_list));
title(str);
str = sprintf('hist_visitage_%s', output_prefix);
utility_save_high_quality_fig(gcf, output_folder, str);

num_total_v = size(age_list, 2);
for age = 30:10:100    
    num_v = size(find(age_list < age), 2);
    fprintf(fp_log, 'num of visit < visiting age %d: %d (%.3f)\n', age, num_v, double(num_v)/double(num_total_v));    
end
cdfplot(age_list);
str = sprintf('cdfplot_visitage_%s', output_prefix);
utility_save_high_quality_fig(gcf, output_folder, str);

%% baseline_age
figure,
hist(baseline_age_list, [0:10:100]);
str = sprintf('histogram of baseline ages (mean:%.1f, std=%.1f)', mean(baseline_age_list), std(baseline_age_list));
title(str);
str = sprintf('hist_baseline_age_%s', output_prefix);
utility_save_high_quality_fig(gcf, output_folder, str);

cdfplot(baseline_age_list);
str = sprintf('cdfplot_baseline_age_%s', output_prefix);
utility_save_high_quality_fig(gcf, output_folder, str);
