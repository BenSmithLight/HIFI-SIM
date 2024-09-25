function HIFI_from_Multiple_single_images(root_path, image_num)
    %% 读取TIF文件
    % clear, clc, close all
    % addpath(genpath('./XxUtils'));
    % addpath(genpath('./Main_fun'));

    % 文件参数
    % filename = '../rawdata/raw_left_2/1.tif';
    % tifname = [filename(1:end - 4), '.tif'];
    filename = [root_path, num2str(image_num), '.tif'];
    outputname = strrep(filename, 'rawdata', 'processed_data');
    % tifname = [outputname(1:end - 4), '.tif'];
    outputname = [outputname(1:end - 4), '_HiFi.tif'];
    tifname = [outputname(1:end - 4), '_merge.tif'];

    % 如果outputname文件夹不存在，则创建
    if ~exist(outputname(1:end - 4), 'dir')
        mkdir(outputname(1:end - 4));
    end

    % 读取从image_num开始的连续9张图片，并将其合成为一个三维数组
    data = imread(filename);
    for i = 1:8
        data = cat(3, data, imread([root_path, num2str(image_num + i), '.tif']));
    end

    % [header, data] = XxReadMRC(filename);
    % [header, data] = read_mrc(filename);

    % data = reshape(data, header(1), header(2), header(3));

    % 将data向左旋转90度
    % data = rot90(data, 1);

    imwrite(data(:, :, 1), tifname);

    for i = 2:9
        imwrite(data(:, :, i), tifname, 'WriteMode', 'append');
    end

    %% 初始化参数
    param = struct();

    % 基本参数

    % 2D SIM参数
    param.phaOff = 0;
    param.fac = ones(1, 2);
    param.nrBands = 2;
    % 3D SIM参数
    % param.phaOff=0;
    % param.fac=ones(1,3);
    % param.nrBands=3;

    % 其他参数
    param.nrDirs = 3; % 方向数，可以根据实际情况修改
    param.nrPhases = 3; % 相位数，可以根据实际情况修改
    param.micronsPerPixel = 65;  % TODO：关键参数
    param.NA = 1.49;  % TODO：关键参数 CCPs,F_actin:1.41  ER,MTs:1.35
    param.lambda = 525;  % TODO：关键参数
    param.attStrength = 0.9;  % default:0.9
    param.a = 0.6; %  阻尼因子1damping2 factor：β default:1
    w1 = 0.9; % Initial optimization Wiener constant：[0.9-2.5]  default:1.2
    ApoFWHM = 99; % 如果设置为99，则使用默认；否则使用设定值  default:99
    WF_flag = 0;
    SR_flag = 0;

    %% 读取tif文件
    N = param.nrDirs * param.nrPhases;
    info = imfinfo(tifname);

    if numel(info) == N
        NPixel = max(info(1).Width, info(1).Height);
        Iraw0 = zeros(info(1).Height, info(1).Width, N);
        Iraw = zeros(NPixel, NPixel, N);

        for j = 1:N
            Iraw0(:, :, j) = double(imread(tifname, j));

            if info(1).Width == info(1).Height || info(1).Width > info(1).Height
                Iraw(1:info(1).Height, :, j) = Iraw0(1:info(1).Height, :, j);
            else
                Iraw(:, 1:info(1).Width, j) = Iraw0(:, 1:info(1).Width, j);
            end

        end

        % 设置输出参数
        param.Format = info.Format;
        param.filename = tifname;
        param.imgSize = NPixel;
        param.Size1 = info(1).Height;
        param.Size2 = info(1).Width;
        param.Iraw = Iraw;
        param.FlagParameter = 0;

        % 显示加载成功信息
        disp('Raw data loaded successfully.');
    else
        warndlg('Raw data frames is not correct!', '!!warn!!', 'modal');
        return;
    end

    %% 开始重建

    FlagParameter = param.FlagParameter;

    if FlagParameter == 0
        %% Generate approximate OTF/PSF
        Iraw = param.Iraw;
        NPixel = size(Iraw, 1);
        param.micronsPerPixel = param.micronsPerPixel * 10 ^ (-3); % 单位：微米
        param.cyclesPerMicron = 1 / (NPixel * param.micronsPerPixel);

        param.attStrength = 0;
        param.OtfProvider = SimOtfProvider(param, param.NA, param.lambda, 1);

        psf = abs(otf2psf((param.OtfProvider.otf)));

        %% Preprocessing
        N = param.nrDirs * param.nrPhases;
        Temp = importImages(Iraw);
        IIraw = deconvlucy(Temp, psf, 5);

        for I = 1:N
            IIrawFFT(:, :, I) = FFT2D(IIraw(:, :, I), false);
        end

        WF = zeros(NPixel, NPixel);
        Tdeconv = zeros(NPixel, NPixel);
        WFdeconv = zeros(NPixel, NPixel, param.nrDirs);
        WFdeconvFFT = zeros(NPixel, NPixel, param.nrDirs);

        for i = 1:param.nrDirs

            for j = 1:param.nrPhases
                Tdeconv = Tdeconv + IIraw(:, :, (i - 1) * param.nrDirs + j);
                WF(:, :) = WF(:, :) + Iraw(:, :, (i - 1) * param.nrDirs + j);
            end

            WFdeconv(:, :, i) = Tdeconv / param.nrPhases;
            WFdeconvFFT(:, :, i) = FFT2D(WFdeconv(:, :, i), false);
        end

        WF = WF / N;
        WF = importImages(WF);

        %% Wide-field
        Size1 = 2 * param.Size1;
        Size2 = 2 * param.Size2;
        WF2 = zeros(Size1, Size2);
        Temp = zeros(2 * size(WF, 1), 2 * size(WF, 2));
        fftWF = fftshift(fft2(WF));
        Temp(size(WF, 1) / 2 + 1:size(WF, 1) / 2 + size(WF, 1), size(WF, 2) / 2 + 1:size(WF, 2) / 2 + size(WF, 2)) = fftWF;
        Temp = abs(ifft2(Temp));
        Temp = 255 * Temp / max(max(Temp));
        WF2(1:Size1, 1:Size2) = Temp(1:Size1, 1:Size2);
        WF2 = importImages(WF2);

        %% Parameter Estimation

        % Kc region MASK
        cnt = [NPixel / 2 + 1, NPixel / 2 + 1];
        param.cutoff = 1000 / (0.5 * param.lambda / param.NA);
        [x, y] = meshgrid(1:NPixel, 1:NPixel);
        rad = sqrt((y - cnt(1)) .^ 2 + (x - cnt(2)) .^ 2);
        Mask = double(rad <= 1.0 * (param.cutoff / param.cyclesPerMicron + 1));
        NotchFilter0 = getotfAtt(NPixel, param.OtfProvider.cyclesPerMicron, 0.5 * param.cutoff, 0, 0);
        NotchFilter = NotchFilter0 .* Mask;
        Mask2 = double(rad <= 1.10 * (param.cutoff / param.cyclesPerMicron + 1));
        NotchFilter2 = NotchFilter0 .* Mask2;

        CrossCorrelation = zeros(size(Mask2, 1), size(Mask2, 2), param.nrDirs);
        k0 = zeros(1, param.nrDirs);

        for I = 1:param.nrDirs
            lb = 2;

            if param.nrBands == 2
                hb = 2;
                fb = lb;
            elseif param.nrBands == 3
                hb = 4;
                fb = hb;
            end

            param.phaOff = 0;
            param.fac = ones(1, param.nrBands);
            separateII = separateBands(IIrawFFT(:, :, (I - 1) * param.nrPhases + 1:I * param.nrPhases), param.phaOff, param.nrBands, param.fac);

            SeparateII{1, I} = separateII;

            if param.nrBands == 2
                c0 = separateII(:, :, 1);
                c2 = separateII(:, :, lb);

                c0 = (c0 ./ (max(max(abs(c0)))));
                c2 = (c2 ./ (max(max(abs(c2)))));

                c0 = c0 .* NotchFilter;
                c2 = c2 .* NotchFilter;

                c0 = FFT2D(c0, false);
                c2 = FFT2D(c2, false);
                c2 = c2 .* conj(c0);
                c2 = c2 ./ max(max(c2));

                vec = fftshift(FFT2D(c2, true));
            elseif param.nrBands == 3
                c0 = separateII(:, :, 1);
                c3 = separateII(:, :, hb);

                c0 = c0 ./ (max(max(abs(separateII(:, :, 1))))) .* NotchFilter;
                c3 = c3 ./ (max(max(abs(separateII(:, :, hb))))) .* NotchFilter;

                c0 = FFT2D(c0, false);
                c3 = FFT2D(c3, false);
                c3 = c3 .* conj(c0);
                c3 = c3 ./ max(max(c3));

                vec = fftshift(FFT2D(c3, true));
            end

            CrossCorrelation(:, :, I) = vec;
            %%
            %         temp=vec.*NotchFilter;
            temp = vec .* NotchFilter2;
            temp = log(1 + abs(temp));
            temp = temp ./ max(max(temp));
            %         MIJ.createImage(temp);

            [yPos, xPos] = find(temp == max(max(temp)));
            peak.xPos = xPos(1);
            peak.yPos = yPos(1);
            k0(I) = sqrt((peak.xPos - cnt(1)) ^ 2 + (peak.yPos - cnt(2)) ^ 2);
        end

        Flag = 0;

        if param.nrDirs > 2 % For very few special cases

            if max(k0) - min(k0) > 8
                Flag = 1;
                Kobject = min(k0);
                %             Kobject=208;
                Mask1 = rad >= (Kobject + 1);
                Mask2 = rad <= (Kobject - 1);
            end

        end

        for I = 1:param.nrDirs
            vec = CrossCorrelation(:, :, I);

            if Flag == 1
                vec(Mask1) = 0;
                vec(Mask2) = 0;
            end

            temp = vec .* NotchFilter2;
            temp = log(1 + abs(temp));
            temp = temp ./ max(max(temp));
            %         MIJ.createImage(temp);

            [yPos, xPos] = find(temp == max(max(temp)));
            peak.xPos = xPos(1);
            peak.yPos = yPos(1);

            cntrl = zeros(10, 30);
            overlap = 0.15;
            step = 2.5;
            bn1 = (param.nrBands - 1) * 2;
            kx = (peak.xPos - cnt(2));
            ky = (peak.yPos - cnt(1));

            separateII = SeparateII{1, I};
            [peak, cntrl] = fitPeak(separateII(:, :, 1) / (max(max(abs(separateII(:, :, 1))))), separateII(:, :, fb) / (max(max(abs(separateII(:, :, fb))))), 1, bn1, param.OtfProvider, -kx, -ky, overlap, step, cntrl);

            if lb ~= hb

                if param.nrBands == 2
                    peak.kx = peak.kx * 2;
                    peak.ky = peak.ky * 2;
                end

                p1 = getPeak(separateII(:, :, 1), separateII(:, :, lb), 0, 1, param.OtfProvider, peak.kx / 2, peak.ky / 2, overlap);
                p2 = getPeak(separateII(:, :, 1) / (max(max(abs(separateII(:, :, 1))))), separateII(:, :, hb) / (max(max(abs(separateII(:, :, 1))))), 0, 2, param.OtfProvider, peak.kx, peak.ky, overlap);

                param.Dir(I).px = -peak.kx / 2;
                param.Dir(I).py = -peak.ky / 2;
                param.Dir(I).phaOff = -phase(p1);
                Temp_m1 = abs(p1);
                Temp_m2 = abs(p2);

                Temp_m1(Temp_m1 > 1.0) = 1;
                Temp_m2(Temp_m2 > 1.0) = 1.0;
                param.Dir(I).modul(1) = Temp_m1;
                param.Dir(I).modul(2) = Temp_m2;
            end

            if lb == hb
                p1 = getPeak(separateII(:, :, 1), separateII(:, :, lb), 1, lb, param.OtfProvider, peak.kx, peak.ky, overlap);
                param.Dir(I).px = -peak.kx;
                param.Dir(I).py = -peak.ky;
                param.Dir(I).phaOff = -phase(p1);
                Temp_m = abs(p1);
                Temp_m(Temp_m > 1.0) = 1.0;
                param.Dir(I).modul = Temp_m;
            end

            K0(I) = sqrt((param.Dir(I).px) ^ 2 + (param.Dir(I).py) ^ 2);
            %%
            %             fittedPeak=vec;
            %             for x=1:30;
            %                 for y=1:10
            %                     fittedPeak(y,x)=cntrl(y,x);
            %                 end
            %             end

            %             figure,imshow(temp,[]);
            %             colormap('hot');
            %             hold on;
            %
            %             if lb~=hb
            %                 f=(param.nrBands-1)/2;
            %             else
            %                 f=1;
            %             end
            %             t=0:0.1:2.1*pi;
            %             x=cnt(2)-peak.kx+10*sin(t);
            %             y=cnt(1)-peak.ky+10*cos(t);
            %             plot(x,y,'-w','LineWidth',2);
            % %             title(['Orientation',num2str(I)]);
            %             %         plotbrowser('on');
            % %         end
            %         %     plotbrowser('off');
        end

        %%
        SIMparam = zeros(3, 5);

        if param.nrPhases == 3

            for i = 1:param.nrDirs
                SIMparam(i, 1) = atan(param.Dir(i).py / param.Dir(i).px) * 180 / pi;
                SIMparam(i, 2) = sqrt((param.Dir(i).px) ^ 2 + (param.Dir(i).py) ^ 2);
                SIMparam(i, 3) = param.Dir(i).phaOff;
                SIMparam(i, 4) = param.Dir(i).modul;

                if param.Dir(i).modul < 0.35
                    param.Dir(i).modul = 0.7;
                end

            end

        elseif param.nrPhases == 5

            for i = 1:param.nrDirs
                SIMparam(i, 1) = atan(param.Dir(i).py / param.Dir(i).px) * 180 / pi;
                SIMparam(i, 2) = sqrt((param.Dir(i).px) ^ 2 + (param.Dir(i).py) ^ 2) * 2;
                SIMparam(i, 3) = param.Dir(i).phaOff;
                SIMparam(i, 4) = param.Dir(i).modul(1);
                SIMparam(i, 5) = param.Dir(i).modul(2);

                if param.Dir(i).modul(1) < 0.35
                    param.Dir(i).modul(1) = 0.7;
                end

                if param.Dir(i).modul(2) < 0.35
                    param.Dir(i).modul(2) = 0.7;
                end

            end

        end

    end

    %% Reconstruction
    siz = size(Iraw(:, :, 1));
    w = siz(2);
    h = siz(1);

    param.attFWHM = 1.0;
    param.OtfProvider = SimOtfProvider(param, param.NA, param.lambda, param.a);

    %% HiFi-SIM：Spectrum optimization
    fftDirectlyCombined = zeros(h * 2, w * 2);

    for I = 1:param.nrDirs
        par = param.Dir(I);
        param.fac(2:param.nrBands) = param.Dir(I).modul(1:param.nrBands - 1);
        param.fac(2:param.nrBands) = param.Dir(I).modul(1:param.nrBands - 1);
        separate = separateBands(IIrawFFT(:, :, (I - 1) * param.nrPhases + 1:I * param.nrPhases), par.phaOff, param.nrBands, param.fac);

        shifted = zeros(2 * h, 2 * w, param.nrPhases);
        shifted(:, :, 1) = placeFreq(separate(:, :, 1));

        for b = 2:param.nrBands
            pos = b * 2 - 2;
            neg = b * 2 - 1;
            shifted(:, :, pos) = placeFreq(separate(:, :, pos));
            shifted(:, :, neg) = placeFreq(separate(:, :, neg));

            shifted(:, :, pos) = NfourierShift(shifted(:, :, pos), - (b - 1) * par.px, - (b - 1) * par.py);
            shifted(:, :, neg) = NfourierShift(shifted(:, :, neg), (b - 1) * par.px, (b - 1) * par.py);
        end

        shifted(:, :, 1) = applyOtf(shifted(:, :, 1), param.OtfProvider, 1, 0, 0, 1, 0);

        for b = 2:param.nrBands
            pos = b * 2 - 2;
            neg = b * 2 - 1;
            shifted(:, :, pos) = applyOtf(shifted(:, :, pos), param.OtfProvider, b, - (b - 1) * par.px, - (b - 1) * par.py, 1, 0);
            shifted(:, :, neg) = applyOtf(shifted(:, :, neg), param.OtfProvider, b, (b - 1) * par.px, (b - 1) * par.py, 1, 0);
        end

        for J = 1:param.nrBands * 2 - 1
            fftDirectlyCombined = fftDirectlyCombined + shifted(:, :, J);
        end

    end

    % Temp1=real(ifft2(fftshift((fftDirectlyCombined))));
    % Temp1(Temp1<0)=0;
    % MIJ.createImage(Temp1);

    w2 = 0.1;
    %%
    param.cutoff = 1000 / (0.5 * param.lambda / param.NA);
    param.sampleLateral = ceil(param.cutoff / param.cyclesPerMicron) + 1;
    K = max([ceil(K0)]);

    if param.nrBands == 2
        cutoff = floor(1 * K) / param.sampleLateral + 1.0;
        R = K;
    elseif param.nrBands == 3
        cutoff = floor(2 * K) / param.sampleLateral + 1.0;
        R = 2 * K;
    end

    otfHiFi = zeros(2 * h, 2 * w);
    otfHiFi = writeApoVector(otfHiFi, param.OtfProvider, cutoff); % Ideal OTF
    Mask = zeros(2 * h, 2 * w);
    Mask(otfHiFi ~= 0) = 1;

    %% Traditional Wiener-SIM
    % if size(Iraw,3)==9
    %     wFilter0=WienerFilterWiener_3D(param);
    % else
    %     wFilter0=WienerFilterWiener_3D(param);
    % end
    % Wk0=otfHiFi./(wFilter0.wDenom+w2^2);
    % fftWiener=real(ifft2(fftshift((fftDirectlyCombined.*Wk0.*Mask))));
    % Temp=fftWiener;
    % Temp(Temp<0)=0;
    % Wiener=255*Temp/max(max(Temp));
    % MIJ.createImage(Wiener);                  % The reconstruction results are displayed in the imageJ window
    % % figure, imshow(Wiener,[]);    % The reconstruction results are displayed in the matlab window
    % % colormap('hot');

    %% HiFi-SIM
    % Step 1
    % if size(Iraw,3)==9
    if param.nrBands == 2
        wFilter1 = WienerFilterW1_2D(param);
    else
        wFilter1 = WienerFilterW1_3D(param);
    end

    Wk1 = otfHiFi ./ (wFilter1.wDenom + w1 ^ 2);

    fftInitialHiFi = fftDirectlyCombined .* Wk1 .* Mask;
    % Temp3=real(ifft2(fftshift((fftInitialHiFi))));
    % Temp3(Temp3<0)=0;
    % MIJ.createImage(Temp3);

    % Step 2
    if size(Iraw, 3) == 9
        wFilter2 = WienerFilterW2_2D(param);
    else
        wFilter2 = WienerFilterW2_3D(param);
    end

    %
    if ApoFWHM == 99
        ApoFWHM = 0.5 * (cutoff - 1);
        ApoFWHM = min(0.5, round(ApoFWHM * 100) / 100);
    end

    apo = apodize_gauss([2 * h, 2 * w], struct('rad', ApoFWHM));
    Wk2 = apo ./ (wFilter2.wDenom + w2 ^ 2);
    fftHiFi = real(ifft2(fftshift((fftInitialHiFi .* Wk2 .* Mask))));

    % Results of wide field
    if WF_flag == 1
        WF = WF2;

        if size(WF, 1) ~= size(WF, 2)
            WF = importImages2(WF);
        end

        %     MIJ.createImage(WF); % The reconstruction results are displayed in the imageJ window
        figure, imshow(WF, []); % The reconstruction results are displayed in the matlab window
        %     colormap('hot');
        param.WF2 = WF;
    end

    % Results of HiFi-SIM
    Size1 = 2 * param.Size1;
    Size2 = 2 * param.Size2;
    HiFi = zeros(Size1, Size2);
    Temp = fftHiFi;
    Temp(Temp < 0) = 0;
    Temp = 255 * Temp / max(max(Temp));
    HiFi(1:Size1, 1:Size2) = Temp(1:Size1, 1:Size2);
    HiFi = importImages2(HiFi);

    if SR_flag == 1
        %     MIJ.createImage(HiFi); % The reconstruction results are displayed in the imageJ window
        figure, imshow(HiFi, []); % The reconstruction results are displayed in the matlab window
        %     colormap('hot');
    end

    param.HiFi = HiFi;

    %% 保存结果
    HiFi_save = uint16(HiFi);
    imwrite(HiFi_save, outputname);
end
