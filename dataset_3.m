clc
clear

% 分块大小
bs = 3;

% 嵌入数据生成
num_D = 1000000;
rng(0);
D = round(rand(1,num_D)*1);

b = 10;

I = imread('./测试图像/Lena.tiff');
origin_I = double(I);
imwrite(uint8(origin_I),'origin.png');

inter_I = origin_I(1:510,1:510);
[row,col] = size(inter_I);

double_I = [inter_I, inter_I];

[encrypt_I, emd_D] = encrypt(double_I, bs, b, D, row, col);
imwrite(uint8(encrypt_I),'encrypt.png');

[recover_I, ext_D]= new_decrypt(encrypt_I, bs, b, row, col);
imwrite(uint8(recover_I),'decrypt.png');

num_emd = length(emd_D);
[MSE, PSNR, ER, BPP] = caulate(row, col, inter_I, recover_I, num_emd);


disp(['嵌入比特数: ' num2str(num_emd) ' bits'] )
disp(['嵌入率: ' num2str(BPP) ' bpp'])
disp(['错误率: ' num2str(ER*100) ' %'])

