function func_plot_datasheet_statistics(output_folder, output_prefix, num_visit_list, trace_year_list, data1_list, data2_list, baseline_age_list, age_list)

global fp_log;


str = sprintf('\nVisit Count: mean:%.1f, std=%.1f\n', mean(num_visit_list), std(num_visit_list));
fprintf(fp_log, str);
figure,
max_value = max(num_visit_list(:));
hist(num_visit_list, [1:1:max_value]);
str = sprintf('Histogram of Visit Count (mean:%.1f, std=%.1f)', mean(num_visit_list), std(num_visit_list));
title(str);
str = sprintf('%s/hist_visit_count_%s', output_folder, output_prefix);
saveas(gcf, str, 'png');
close(gcf);

str = sprintf('Tracing Year: mean:%.1f, std=%.1f\n', mean(trace_year_list), std(trace_year_list));
fprintf(fp_log, str);
figure,
max_value = max(trace_year_list(:));
hist(trace_year_list, [1:1:max_value]);
str = sprintf('Histogram of Tracing Years (mean:%.1f, std=%.1f)', mean(trace_year_list), std(trace_year_list));
title(str);
str = sprintf('%s/hist_tracing_year_%s', output_folder, output_prefix);
saveas(gcf, str, 'png');
close(gcf);

str = sprintf('data1: mean:%.1f, std=%.1f\n',mean(data1_list), std(data1_list));
fprintf(fp_log, str);
figure,
%hist(VFI_list, [0:1:100]);
num_bins = 20;
hist(data1_list, num_bins);

str = sprintf('Histogram of data1 (mean:%.1f, std=%.1f)', mean(data1_list), std(data1_list));
title(str);
%xlim([0 100]);
str = sprintf('%s/hist_data1_%s', output_folder, output_prefix);
saveas(gcf, str, 'png');
close(gcf);

num_data1_data = size(data1_list, 1);
%fprintf(fp_log, 'num VFI data = %d\n', num_VFI_data);

%for VFI_value = 100:(-1):85
%    fprintf(fp_log, 'VFI == %d (%.2f %%)\n', VFI_value, size(find(VFI_list == VFI_value), 1)/double(num_VFI_data) * 100.0);    
%end
%fprintf(fp_log, '\n');
%for VFI_value = 95:(-5):5
    %fprintf(fp_log, 'VFI >=  %d (%.2f %%)\n', VFI_value, size(find(VFI_list >= VFI_value), 1)/double(num_VFI_data) * 100.0);
%end

%==================================

str = sprintf('\ndata2: mean:%.1f, std=%.1f\n', mean(data2_list), std(data2_list));
fprintf(fp_log, str);
figure,
%hist(RNFL_list, [0:5:155]);
num_bins = 20;
hist(data2_list, num_bins);

str = sprintf('Histogram of data2 (mean:%.1f, std=%.1f)', mean(data2_list), std(data2_list));
title(str);
%xlim([0 155]);
fprintf(fp_log, 'min data2 (%.1f)\n', min(data2_list));
fprintf(fp_log, 'max data2 (%.1f)\n', max(data2_list));
str = sprintf('%s/hist_data2_%s', output_folder, output_prefix);
saveas(gcf, str, 'png');
close(gcf);

str = sprintf('\nVisit ages: mean:%.1f, std=%.1f\n', mean(age_list), std(age_list));
fprintf(fp_log, str);
figure,
hist(age_list, [0:10:110]);
str = sprintf('Histogram of Visit Ages (mean:%.1f, std=%.1f)', mean(age_list), std(age_list));
title(str);
xlabel('age');
str = sprintf('%s/hist_visitage_%s', output_folder, output_prefix);
saveas(gcf, str, 'png');
close(gcf);

num_total_v = length(age_list);
for age = 30:10:100    
    num_v = size(find(age_list <= age), 1);
    fprintf(fp_log, 'num of visit <= visiting age %d: %d (%.1f%%)\n', age, num_v, double(num_v)/double(num_total_v) * 100);    
end
figure,
cdfplot(age_list);
str = sprintf('%s/cdfplot_visitage_%s', output_folder, output_prefix);
saveas(gcf, str, 'png');
close(gcf);

%% baseline_age
str = sprintf('baseline age: mean:%.1f, std=%.1f\n', mean(baseline_age_list), std(baseline_age_list));
fprintf(fp_log, str);
figure,
hist(baseline_age_list, [0:10:110]);
str = sprintf('Histogram of Baseline Ages (mean:%.1f, std=%.1f)', mean(baseline_age_list), std(baseline_age_list));
title(str);
xlabel('age');
str = sprintf('%s/hist_baseline_age_%s', output_folder, output_prefix);
saveas(gcf, str, 'png');
close(gcf);

figure,
cdfplot(baseline_age_list);
str = sprintf('%s/cdfplot_baseline_age_%s', output_folder, output_prefix);
saveas(gcf, str, 'png');
close(gcf);
