clc
clear

% 分块大小
bs = 3;

% 嵌入数据生成
num_D = 1000000;
rng(0);
D = round(rand(1,num_D)*1);

% 每块嵌入5 bit数据
b = 5;

I = imread('./测试图像/Lena.tiff');
origin_I = double(I);

inter_I = origin_I(1:510,1:510);
[pre_I,num_Of,Overflow] = predict_error(inter_I);

[row,col] = size(inter_I);

double_I = [pre_I, pre_I];

figure(1);
imshow(double_I,[]);title('拼接图像');
%% 加密 & 解密

[encrypt_I, emd_D] = encrypt(double_I, bs, b, D, row, col);

[recover_I, ext_D]= decrypt(encrypt_I, bs, b, row, col);

num_emd = length(emd_D);
[MSE, PSNR, ER, BPP] = caulate(row, col, inter_I, recover_I, num_emd);

figure(1);
subplot(221);
imshow(inter_I,[]);title('原始图像');
subplot(222);
imshow(double_I,[]);title('拼接图像');
subplot(223);
imshow(encrypt_I,[]);title('加密图像');
subplot(224);
imshow(recover_I,[]);title('复原图像');




% disp(['嵌入比特数: ' num2str(num_emd) ' bits'] )
% disp(['嵌入率: ' num2str(BPP) ' bpp'])
% disp(['错误率: ' num2str(ER) ])


