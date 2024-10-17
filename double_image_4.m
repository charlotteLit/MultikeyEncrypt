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


[MSE, PSNR, ER, BPP] = caulate(row, col, inter_I, recover_I);

check = isequal(emd_D,ext_D);
if check == 1  
    disp('提取数据与嵌入数据完全相同！')
    disp(['嵌入比特数: ' num2str(num_emD) ' bits'] )
    disp(['嵌入率: ' num2str(bpp) ' bpp'])
else
    disp('Warning！数据提取错误！')
end