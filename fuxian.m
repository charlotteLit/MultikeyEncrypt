clc
clear

% 分块大小
bs = 5;

% 嵌入数据生成
num_D = 1000000;
rng(0);
D = round(rand(1,num_D)*1);

% 每块嵌入5 bit数据
b = 5;


%% 图像数据集信息
I_file_path = '.\测试图像\'; %测试图像数据集文件夹路径
I_path_list = dir(strcat(I_file_path,'*.tiff'));  %获取该文件夹中所有pgm格式的图像
img_num = length(I_path_list); %获取图像总数量

%% 运行信息

for i=1:img_num
    %-------------------------------读取图像-------------------------------%
    I_name = I_path_list(i).name; %图像名
    I = imread(strcat(I_file_path,I_name));%读取图像
    origin_I = double(I);

    inter_I = origin_I(1:510,1:510);
    % inter_I = origin_I;

    [row,col] = size(inter_I);

    [encrypt_I, emd_D] = fuxian_encrypt(inter_I, bs, b, D, row, col);

    [recover_I, ext_D]= fuxian_decrypt(encrypt_I, bs, b, row, col);

    num_emd = length(emd_D);
    [MSE, PSNR, ER, BPP] = caulate(row, col, inter_I, recover_I, num_emd);

    disp(['image:' I_name])
    disp(['嵌入比特数: ' num2str(num_emd) ' bits'] )
    disp(['嵌入率: ' num2str(BPP) ' bpp'])
    disp(['错误率: ' num2str(ER*100) '%'])

end