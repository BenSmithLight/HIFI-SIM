filename1 = '../../Picture/Pre_process/ER/RawGTSIMData/RawGTSIMData_level_01.mrc';
filename2 = '../../Picture/Pre_process/ER/RawGTSIMData/RawGTSIMData_level_02.mrc';

[header1, data1] = read_mrc(filename1);

data1 = reshape(data1, header1(1), header1(2), header1(3));

[header2, data2] = read_mrc(filename2);

data2 = reshape(data2, header2(1), header2(2), header2(3));

data1_1 = data1(:,:,1);
data2_1 = data2(:,:,1);

figure;  imshow(data1(:,:,1), []); title('Level 1');
figure;  imshow(data2(:,:,1), []); title('Level 2');
