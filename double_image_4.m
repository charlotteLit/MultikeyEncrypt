clc
clear

% �ֿ��С
bs = 4;

% Ƕ����������
num_D = 1000000;
rng(0);
D = round(rand(1,num_D)*1);

% ÿ��Ƕ��5 bit����
b = 5;

I = imread('./����ͼ��/Lena.tiff');
origin_I = double(I);

inter_I = origin_I;
[row,col] = size(inter_I);

double_I = [inter_I, inter_I];

[encrypt_I, emd_D] = encrypt(double_I, bs, b, D, row, col);

[recover_I, ext_D]= decrypt(encrypt_I, bs, b, row, col);

num_emd = length(emd_D);
[MSE, PSNR, ER, BPP] = caulate(row, col, inter_I, recover_I, num_emd);


disp(['Ƕ�������: ' num2str(num_emd) ' bits'] )
disp(['Ƕ����: ' num2str(BPP) ' bpp'])
disp(['������: ' num2str(ER) ' bpp'])