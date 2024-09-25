close all;clear;clc;
currPath = fileparts(mfilename('fullpath'));    
cd(currPath);                                   
addpath(genpath('./Main_fun'));
% addpath(genpath('./XxUtils'));

% raw_left_2
root_path = '../rawdata/raw_right_2/';
image_num = 1;
filename = [root_path, num2str(image_num), '.tif'];

% while exist(filename, 'file')
%     disp(['Processing image ', num2str(image_num)]);
%     main_BioSR(filename);
%     image_num = image_num + 9;
% end

% 在上面的基础上，要判断所有9张图片是否都存在，如果有一张不存在，就不处理
while exist(filename, 'file')
    disp(['Processing image ', num2str(image_num)]);
    for i = 0:8
        filename = [root_path, num2str(image_num + i), '.tif'];
        if ~exist(filename, 'file')
            disp(['Image ', num2str(image_num + i), ' does not exist!']);
            break;
        end
    end
    if i == 8
        HIFI_from_Multiple_single_images(root_path, image_num);
    end
    image_num = image_num + 9;
end