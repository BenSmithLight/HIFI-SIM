clear; clc; close all;

filename = '../../Picture/Pre_process/ER/RawSIMData/RawSIMData_level_01.mrc';

for i = 1:6
    if i < 10
        filename_process = strrep(filename, '01', ['0', num2str(i)]);
    else
        filename_process = strrep(filename, '01', num2str(i));
    end
    disp(['Processing: ', filename_process]);
    outputname = strrep(filename_process, 'Pre', 'After');
    tifname = [outputname(1:end - 4), '.tif'];
    outputname = [outputname(1:end - 4), '_HiFi.tif'];

    [header, data] = XxReadMRC(filename_process);

    data = reshape(data, header(1), header(2), header(3));

    imwrite(data(:, :, 1), tifname);

    for j = 2:9
        imwrite(data(:, :, i), tifname, 'WriteMode', 'append');
    end
end
