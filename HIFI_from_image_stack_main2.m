close all;clear;clc;
currPath = fileparts(mfilename('fullpath'));
cd(currPath);
addpath(genpath('./Main_fun'));
% addpath(genpath('./XxUtils'));

% full size path
% filename = '../rawdata/lifeact-exocytosis/EGFP/5_488_Em525_TIRFSIM-2_1.7-INS1-Lifeact-EGFP-Vamp2-pHmScarlet-20240415.tif';
% filename = '../rawdata/lifeact-exocytosis/mScarlet/10_561_Em609_TIRFSIM-2_1.7-INS1-Lifeact-mScarlet-Vamp2-pHlourin-20240415.tif';

% SN2N path
% filename = '../rawdata/lifeact_SN2N/raw1_left/raw1_left.tif';
filename = '../rawdata/lifeact_SN2N/raw1_right/raw1_right.tif';
% filename = '../rawdata/lifeact_SN2N/raw2_left/raw2_left.tif';
% filename = '../rawdata/lifeact_SN2N/raw2_right/raw2_right.tif';

image_frame = 1;

while image_frame < 1800
    disp(['Processing image ', num2str(image_frame), '-', num2str(image_frame + 8), '...']);

    origin_test(filename, image_frame);

    image_frame = image_frame + 9;
end
