% ------------------------------------------------------------------------
% 该函数用于生成DN模块的训练数据
% ------------------------------------------------------------------------

function [data_seg, gt_seg] = XxDataSeg_ForTrain_rDL(data, gt, num_seg, Nx, Ny, RotFlag)

% ------------------------------------------------------------------------
% Randomly crop and rotate image patch pairs from data and gt 
% 
% usage:  [DATA_SEG, GT_SEG] = XxDataSeg_ForTrain(data, gt, num_seg, Nx, 
%                              Ny, RotFlag)
% where,
%    DATA_SEG    -- image patches of raw SIM data (input of networks)
%    GT_SEG      -- image patches of SIM SR data (ground truth)
%    data        -- raw SIM data before segmentation
%    gt          -- SIM SR data before segmentation
%    num_seg     -- number of patches to crop
%    Nx, Ny      -- patch size to crop
%    RotFlag     -- >=1 for random angle rotation, 0 for no rotation
%
% Author: Chang Qiao
% Email: qc17@mails.tsinghua.edu.cn
% Version  : 2020/5/15
% ------------------------------------------------------------------------

addpath(genpath('./XxUtils'));

[r,l,~] = size(data);
[r_h,~,~] = size(gt);

if RotFlag == 0
    new_r = Nx;
    new_l = Ny;
else
    new_r = ceil(Nx * 1.5); % ceil即向上取整
    new_l = ceil(Ny * 1.5);
    if new_r > r || new_l > l % 如果patch尺寸大于图像尺寸，则不旋转
        new_r = r;
        new_l = l;
        RotFlag = 0;
    end
end

% calculate foreground mask 计算前景掩膜
thresh_mask = 1e-2;
ksize = 3; 
% gt_cal_mask = imresize(gt, 1/sr_ratio, 'bicubic'); % 由于gt改为三维，下面对第三维度求和，所以这里不需要resize
gt_cal_mask = sum(gt,3); % 对gt的第三维度求和，得到二维图像
mask = XxCalMask(gt_cal_mask,ksize,thresh_mask); % 计算前景掩膜，要求输入gt_cal_mask为二维图像
% 如果前景掩膜像素点数小于1000，则降低阈值
while sum(mask(:)) < 1e3
    thresh_mask = thresh_mask * 0.8;
    mask = XxCalMask(gt_cal_mask,ksize,thresh_mask);
end

% figure(1);
% subplot(1,2,1), imshow(gt,[]);
% subplot(1,2,2), imshow(mask,[]);

Y = 1:l;
X = 1:r;
[X,Y] = meshgrid(Y,X); % 生成网格点坐标矩阵
point_list = zeros(sum(mask(:)),2); % 生成前景掩膜中的点坐标
point_list(:,1) = Y(mask(:)); % mask(:)将mask矩阵转换为列向量
point_list(:,2) = X(mask(:)); % point_list为前景掩膜中的点坐标
l_list = size(point_list,1); % 前景掩膜中的点数

halfx = round(new_r / 2); % 1.5倍patch的一半
halfy = round(new_l / 2);

data_seg = [];
gt_seg = [];

%% Crop patches 裁剪patch
thresh_ar = 0;
thresh_ar_gt = 0;
thresh_sum_gt = 0;
count = 0;

for i = 1:num_seg
    p = randi(l_list,1); % 随机选取一个前景掩膜中的点
    y1 = point_list(p, 1) - halfy + 1; % patch左上角坐标
    y2 = point_list(p, 1) + halfy; % patch右下角坐标
    x1 = point_list(p, 2) - halfx + 1; % patch左上角坐标
    x2 = point_list(p, 2) + halfx; % patch右下角坐标
    count_rand = 0; % 随机选取patch的次数
    % 如果patch超出图像范围，则重新随机选取patch
    while (y1<1 || y2>r || x1<1 || x2>l)
        p = randi(l_list,1);
        y1 = point_list(p, 1) - halfy + 1;
        y2 = point_list(p, 1) + halfy;
        x1 = point_list(p, 2) - halfx + 1;
        x2 = point_list(p, 2) + halfx;
        count_rand = count_rand + 1;
        if count_rand > 1e3, break; end % 如果随机选取patch的次数超过1000，则跳出循环
    end
    if count_rand > 1e3, break; end % 如果随机选取patch的次数超过1000，则跳出循环
    
    if RotFlag >= 1 % if random rotate 如果随机旋转
        degree = randi(360, 1); % 获取随机旋转角度

        % 对data中选定区域第三维求和，得到二维图像，然后旋转，使用双线性插值，并让旋转后的图像大小与patch大小相同
        patch = imrotate(sum(data(x1:x2,y1:y2,:), 3),degree,'bilinear','crop'); % 192*192 double
        tx = round(new_r/2)-round(Nx/2); % patch的起始坐标
        ty = round(new_l/2)-round(Ny/2);
        patch = patch(tx+1:tx+Nx,ty+1:ty+Ny); % 128*128
        active_range = double(prctile(patch(:),99.9)) / double(prctile(patch(:),0.1)+1e-2); % 计算图像像素值的99.9百分位数和0.1百分位数的比值，这可以用来评估图像的对比度
        
        % patch_gt = imrotate(gt(sr_ratio*x1-1:sr_ratio*x2,sr_ratio*y1-1:sr_ratio*y2,1),degree,'bilinear','crop'); % 384*384 使用了sr_ratio来调整区域的大小
        % tx_gt = round(new_r*sr_ratio/2)-round(Nx*sr_ratio/2);  % gt_patch的起始坐标
        % ty_gt = round(new_l*sr_ratio/2)-round(Ny*sr_ratio/2);
        % patch_gt = patch_gt(tx_gt+1:tx_gt+sr_ratio*Nx,ty_gt+1:ty_gt+sr_ratio*Ny); % 256*256
        % sum_patch_gt = sum(patch_gt(:));
        % active_range_gt = double(prctile(patch_gt(:),99.9)) / double(prctile(patch_gt(:),0.1)+1); % 对比度评估
        % 由于gt改为三维，因此需要二维输入的，均对第三维求和；并且gt大小与data大小相同，无需缩放
        patch_gt = imrotate(sum(gt(x1:x2,y1:y2,:), 3),degree,'bilinear','crop'); % 192*192
        tx_gt = round(new_r/2)-round(Nx/2); % gt_patch的起始坐标
        ty_gt = round(new_l/2)-round(Ny/2);
        patch_gt = patch_gt(tx_gt+1:tx_gt+Nx,ty_gt+1:ty_gt+Ny); % 128*128
        sum_patch_gt = sum(patch_gt(:));
        active_range_gt = double(prctile(patch_gt(:),99.9)) / double(prctile(patch_gt(:),0.1)+1); % 对比度评估
    else % not random rotate % TODO 这里没有做修改
        patch = sum(data(x1:x2,y1:y2,:), 3);
        active_range = double(prctile(patch(:),99.9)) / double(prctile(patch(:),0.1)+1);
        
        % patch_gt = gt(sr_ratio*x1-1:sr_ratio*x2,sr_ratio*y1-1:sr_ratio*y2,1);
        patch_gt = sum(gt(x1:x2,y1:y2,:), 3);
        sum_patch_gt = sum(patch_gt(:));
        active_range_gt = double(prctile(patch_gt(:),99.9)) / double(prctile(patch_gt(:),0.1)+1);
    end
    
    % 如果gt_patch的对比度评估值小于阈值、gt_patch的像素和小于阈值、patch的对比度评估值小于阈值，则重新随机选取patch
    while active_range_gt < thresh_ar_gt || sum_patch_gt < thresh_sum_gt || active_range < thresh_ar
        x1 = randi(r - new_r + 1, 1);
        x2 = x1 + new_r - 1;
        y1 = randi(l - new_l + 1, 1);
        y2 = y1 + new_l - 1;
        
        if RotFlag >= 1 % if random rotate
            degree = randi(360, 1);
            patch = imrotate(sum(data(x1:x2,y1:y2,:), 3),degree,'bilinear','crop');
            tx = round(new_r/2)-round(Nx/2);
            ty = round(new_l/2)-round(Ny/2);
            patch = patch(tx+1:tx+Nx,ty+1:ty+Ny);
            active_range = double(prctile(patch(:),99.9)) / double(prctile(patch(:),0.1)+1e-2);
            
            % patch_gt = imrotate(gt(sr_ratio*x1-1:sr_ratio*x2,sr_ratio*y1-1:sr_ratio*y2,1),degree,'bilinear','crop');
            % tx_gt = round(new_r*sr_ratio/2)-round(Nx*sr_ratio/2);
            % ty_gt = round(new_l*sr_ratio/2)-round(Ny*sr_ratio/2);
            % patch_gt = patch_gt(tx_gt+1:tx_gt+sr_ratio*Nx,ty_gt+1:ty_gt+sr_ratio*Ny);
            % sum_patch_gt = sum(patch_gt(:));
            % active_range_gt = double(prctile(patch_gt(:),99.9)) / double(prctile(patch_gt(:),0.1)+1);
            patch_gt = imrotate(sum(gt(x1:x2,y1:y2,:), 3),degree,'bilinear','crop'); % 384*384
            tx_gt = round(new_r/2)-round(Nx/2); % gt_patch的起始坐标
            ty_gt = round(new_l/2)-round(Ny/2);
            patch_gt = patch_gt(tx_gt+1:tx_gt+Nx,ty_gt+1:ty_gt+Ny); % 256*256
            sum_patch_gt = sum(patch_gt(:));
            active_range_gt = double(prctile(patch_gt(:),99.9)) / double(prctile(patch_gt(:),0.1)+1); % 对比度评估
        else % not random rotate % TODO 这里没有做修改
            patch = sum(data(x1:x2,y1:y2,:), 3);
            active_range = double(prctile(patch(:),99.9)) / double(prctile(patch(:),0.1)+1e-2);
            
            % patch_gt = gt(sr_ratio*x1-1:sr_ratio*x2,sr_ratio*y1-1:sr_ratio*y2,1);
            patch_gt = sum(gt(x1:x2,y1:y2,:), 3);
            sum_patch_gt = sum(patch_gt(:));
            active_range_gt = double(prctile(patch_gt(:),99.9)) / double(prctile(patch_gt(:),0.1)+1);
        end
        
        count = count + 1;
        % 如果重新随机选取patch的次数超过100，则降低阈值
        if count > 1e2
            thresh_ar_gt = thresh_ar_gt * 0.9;
            thresh_sum_gt = thresh_sum_gt * 0.9;
            thresh_ar = thresh_ar * 0.9;
            count = 0;
        end
        
%         fprintf('ar_gt    = %.2f\n',active_range_gt);
%         fprintf('sum_gt   = %.2f\n',sum_patch_gt);
%         fprintf('ar_data  = %.2f\n',active_range);
%         figure(1);
%         subplot(1,2,1), imagesc(patch(:,:,1)); axis image, axis off;
%         subplot(1,2,2), imagesc(patch_gt(:,:,1)); axis image, axis off;
    end
    
    if RotFlag >= 1
        tdata = imrotate(data(x1:x2,y1:y2,:),degree,'bilinear','crop'); % 192*192*9
        % tgt = imrotate(gt(sr_ratio*x1-1:sr_ratio*x2,sr_ratio*y1-1:sr_ratio*y2,:),degree,'bilinear','crop'); %384*384
        tgt = imrotate(gt(x1:x2,y1:y2,:),degree,'bilinear','crop'); % 192*192*9
        if isempty(data_seg)
            data_seg = tdata(tx+1:tx+Nx,ty+1:ty+Ny,:);  % 128*128*9
        else
            h  = size(tdata, 3);
            data_seg(:,:,end+1:end+h) = tdata(tx+1:tx+Nx,ty+1:ty+Ny,:);
        end
        if isempty(gt_seg)
            % gt_seg = tgt(tx_gt+1:tx_gt+sr_ratio*Nx,ty_gt+1:ty_gt+sr_ratio*Ny,:); % 256*256
            gt_seg = tgt(tx_gt+1:tx_gt+Nx,ty_gt+1:ty_gt+Ny,:); % 128*128*9
        else
            h_gt = size(tgt, 3);
            % gt_seg(:,:,end+1:end+h_gt) = tgt(tx_gt+1:tx_gt+sr_ratio*Nx,ty_gt+1:ty_gt+sr_ratio*Ny,:);
            gt_seg(:,:,end+1:end+h_gt) = tgt(tx_gt+1:tx_gt+Nx,ty_gt+1:ty_gt+Ny,:);
        end
    else  % TODO 这里没有做修改
        if isempty(data_seg)
            data_seg = data(x1:x2,y1:y2,:);
        else
            h  = size(data, 3);
            data_seg(:,:,end+1:end+h) = data(x1:x2,y1:y2,:);
        end
        if isempty(gt_seg)
            % gt_seg = gt(sr_ratio*x1-1:sr_ratio*x2,sr_ratio*y1-1:sr_ratio*y2,:);
            gt_seg = gt(x1:x2,y1:y2,:);
        else
            h_gt = size(gt, 3);
            % gt_seg(:,:,end+1:end+h_gt) = gt(sr_ratio*x1-1:sr_ratio*x2,sr_ratio*y1-1:sr_ratio*y2,:);
            gt_seg(:,:,end+1:end+h_gt) = gt(x1:x2,y1:y2,:);
        end
    end
end

end