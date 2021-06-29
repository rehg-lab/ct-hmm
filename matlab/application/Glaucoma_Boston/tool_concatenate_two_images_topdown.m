%% this tool can combine the two images from the same patient in a top-bottom fashion
%% modify the path for the two input image folder
%% the script reads the image from the in_folder1, and combine the image from the in_folder2 which has the same image name
%% specify the name for your output folder

in_folder1 = 'output_decoding_minvisit5/decoded_traj/';
in_folder2 = 'output_decoding_minvisit5/decoded_traj/';
out_folder = 'combined_img';
mkdir(out_folder);

filestr = sprintf('%s/*.png', folder1);
img_file_list = dir(filestr); 

num_img = size(img_file_list, 1);

for i = 1:num_img

    img_name = img_file_list(i).name;
    
    file1 = sprintf('%s/%s', in_folder1, img_name);
    file2 = sprintf('%s/%s', in_folder2, img_name);

    %file1 = 'output_decoding_minvisit5/decoded_traj/decoded_traj_001_2217_OS.png';
    %file2 = 'output_decoding_minvisit5/decoded_traj/decoded_traj_002_2372_OS.png';

    img1 = imread(file1);
    img2 = imread(file2);

    size_img1 = size(img1);
    size_img2 = size(img2);

    combined_row = size_img1(1) + size_img2(1);
    combined_col = max(size_img1(2), size_img2(2));

    %% create the combined image
    combined_img = uint8(ones(combined_row, combined_col, 3) * 255);

    %% assign image 1 at top
    combined_img(1:size_img1(1), 1:size_img1(2), :) = img1;

    %% assign image 2 at bottom
    begin_row = size_img2(1) + 1;
    end_row = begin_row + size_img2(1) - 1;
    combined_img(begin_row:end_row, 1:size_img2(2), :) = img2;

    %% save the combined image
    out_img_file = sprintf('%s/%s', out_folder, img_name);    
    imwrite(combined_img, out_img_file);

end

