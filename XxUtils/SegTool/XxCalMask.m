function mask = XxCalMask(img, ksize, thresh)

% ------------------------------------------------------------------------
% XxCalMask: calculate the foreground mask of img 计算图像的前景掩码
%
% usage:  mask = XxCalMask(img, ksize, thresh)
% where,
%    img       -- 2D img to calculate foreground 二维图像
%    ksize     -- the size of first gaussian kernel, typically set as 5~10 第一个高斯核的大小，通常设置为5~10
%    thresh    -- lower thresh leads to larger mask, typically set as 
%                 [1e-3, 5e-2] 阈值越低，掩码越大，通常设置为[1e-3, 5e-2]
%
% Author: Chang Qiao
% Email: qc17@mails.tsinghua.edu.cn
% Version: 2020/5/15
% ------------------------------------------------------------------------

% 检查输入参数
if nargin < 3, thresh = 5e-2; end % 如果没有输入阈值，则默认为5e-2
if nargin < 2, ksize = 10; end % 如果没有输入高斯核大小，则默认为10

% 创建一个高斯核，用于图像滤波
kernel = fspecial('gaussian',[ksize,ksize],ksize); % 生成高斯核
fd = imfilter(img,kernel,'replicate'); % 对图像进行高斯滤波，得到前景细节
% 创建一个大的高斯核，用于估计图像背景
kernel = fspecial('gaussian',[100,100],50); % 生成高斯核
bg = imfilter(img,kernel,'replicate'); % 对图像进行高斯滤波，得到背景大致情况

mask = fd - bg; % 计算前景和背景的差异
mask(mask >= thresh) = 1; % 将差异大于等于阈值的部分设置为1，即认为是前景
mask(mask ~= 1) = 0; % 将差异小于阈值的部分设置为0，即认为是背景
mask = logical(mask); % 将mask转换为逻辑类型

end