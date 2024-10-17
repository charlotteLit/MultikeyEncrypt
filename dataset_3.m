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

%% 图像数据集信息(BOSSbase_1.01),格式:PGM,数量:10000；
I_file_path = 'D:\ImageDatabase\BOSSbase_1.01\'; %测试图像数据集文件夹路径
I_path_list = dir(strcat(I_file_path,'*.pgm'));  %获取该文件夹中所有pgm格式的图像
img_num = length(I_path_list); %获取图像总数量
%% 记录每张图像的的相关信息
bpp_BOSSbase = zeros(1,img_num); %记录每张图像的嵌入率
er_BOSSbase = zeros(1,img_num);%记录每张图像的溢出像素个数

for i=1:img_num
    %-------------------------------读取图像-------------------------------%
    I_name = I_path_list(i).name; %图像名
    I = imread(strcat(I_file_path,I_name));%读取图像
    origin_I = double(I);

    inter_I = origin_I(1:510,1:510);
    [row,col] = size(inter_I);

    double_I = [inter_I, inter_I];

    [encrypt_I, emd_D] = encrypt(double_I, bs, b, D, row, col);

    [recover_I, ext_D]= decrypt(encrypt_I, bs, b, row, col);

    [MSE, PSNR, ER, BPP] = caulate(row, col, inter_I, recover_I);
    
    er_BOSSbase(i) = ER;
    bpp_BOSSbase(i) = BPP;

    check = isequal(emd_D,ext_D);
    if check == 1  
        disp('提取数据与嵌入数据完全相同！')
        disp(['嵌入比特数: ' num2str(num_emD) ' bits'] )
        disp(['嵌入率: ' num2str(bpp) ' bpp'])
    else
        disp('Warning！数据提取错误！')
    end
    fprintf(['第 ',num2str(i),' 幅图像-------- OK','\n\n']);
end
