clc
clear

% 分块大小
bs = 4;

% 嵌入数据生成
num_D = 1000000;
rng(0);
D = round(rand(1,num_D)*1);

% 每块嵌入5 bit数据
b = 5;

I = imread('./测试图像/Lena.tiff');
origin_I = double(I);

inter_I = origin_I;
[row,col] = size(inter_I);

double_I = [inter_I, inter_I];

[encrypt_I, emd_D] = encrypt(double_I, bs, b, D, row, col);

[recover_I, ext_D]= decrypt(encrypt_I, bs, b, row, col);

num_emd = length(emd_D);
[MSE, PSNR, ER, BPP] = caulate(row, col, inter_I, recover_I, num_emd);


disp(['嵌入比特数: ' num2str(num_emd) ' bits'] )
disp(['嵌入率: ' num2str(BPP) ' bpp'])
disp(['错误率: ' num2str(ER) ' bpp'])