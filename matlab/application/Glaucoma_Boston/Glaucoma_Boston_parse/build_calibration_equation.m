if (code == (12+1)) % prototype + 1_2
        i = prototype_oct_1_2_count + 1;
        prototype_oct_1_2_data{i}.x = mean_rnfl_prototype;
        prototype_oct_1_2_data{i}.y = mean_rnfl_oct_1_2;
        prototype_oct_1_2_count = i;
    elseif (code == (6+1)) % 1_2 + stratus
        i = oct_1_2_stratus_count + 1;
        oct_1_2_stratus_data{i}.x = mean_rnfl_oct_1_2;
        oct_1_2_stratus_data{i}.y = mean_rnfl_stratus;
        oct_1_2_stratus_count = i;
    elseif (code == (3+1)) % stratus + cirrus
        i = stratus_cirrus_count + 1;
        stratus_cirrus_data{i}.x = mean_rnfl_stratus;
        stratus_cirrus_data{i}.y = mean_rnfl_cirrus;
        stratus_cirrus_count = i;
    end


%%===============================================================================
%% build calibration equations between machines
%% prototype + 1_2
%count = prototype_oct_1_2_count;
%data = prototype_oct_1_2_data;
%outfile = 'regress_prototype_oct_1_2.png';
%[slope, intercept] = func_linear_regress_plot(count, data, outfile);
%prototype_oct_1_2_slope = slope;
%prototype_oct_1_2_intercept = intercept;

str = sprintf('Prototype_Commercial_OCT_Overlap2.xlsx');
[num,txt,raw] = xlsread(str);
num_row = size(num, 1);
num_col = size(num, 2);

%% field we want to use
field_list = {'mean_rnfl_prototype_rnfl', 'mean_rnfl_commercial_rnfl'};              
field_col = zeros(2, 1);

%% find field position
for f = 1:2
    for c = 1:num_col
        if (strcmp(txt(1, c), field_list{f})==1)
            field_col(f) = c;
            break;
        end
    end
end
data = cell(num_row, 1);
for r = 1:num_row
     data{r}.x = num(r, field_col(1));
     data{r}.y = num(r, field_col(2));
end

count = num_row;
outfile = 'regress_prototype_oct_1_2.png';
[slope, intercept] = func_linear_regress_plot(count, data, outfile);
prototype_oct_1_2_slope = slope;
prototype_oct_1_2_intercept = intercept;


%% oct_1_2 + stratus
count = oct_1_2_stratus_count;
data = oct_1_2_stratus_data;
outfile = 'regress_oct_1_2_stratus.png';
[slope, intercept] = func_linear_regress_plot(count, data, outfile);
oct_1_2_stratus_slope = slope;
oct_1_2_stratus_intercept = intercept;

%% stratus + cirrus
count = stratus_cirrus_count;
data = stratus_cirrus_data;
outfile = 'regress_stratus_cirrus.png';
[slope, intercept] = func_linear_regress_plot(count, data, outfile);
stratus_cirrus_slope = slope;
stratus_cirrus_intercept = intercept;
