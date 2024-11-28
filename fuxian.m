clc
clear

% �ֿ��С
bs = 5;

% Ƕ����������
num_D = 1000000;
rng(0);
D = round(rand(1,num_D)*1);

% ÿ��Ƕ��5 bit����
b = 5;


%% ͼ�����ݼ���Ϣ
I_file_path = '.\����ͼ��\'; %����ͼ�����ݼ��ļ���·��
I_path_list = dir(strcat(I_file_path,'*.tiff'));  %��ȡ���ļ���������pgm��ʽ��ͼ��
img_num = length(I_path_list); %��ȡͼ��������

%% ������Ϣ

for i=1:img_num
    %-------------------------------��ȡͼ��-------------------------------%
    I_name = I_path_list(i).name; %ͼ����
    I = imread(strcat(I_file_path,I_name));%��ȡͼ��
    origin_I = double(I);

    inter_I = origin_I(1:510,1:510);
    % inter_I = origin_I;

    [row,col] = size(inter_I);

    [encrypt_I, emd_D] = fuxian_encrypt(inter_I, bs, b, D, row, col);

    [recover_I, ext_D]= fuxian_decrypt(encrypt_I, bs, b, row, col);

    num_emd = length(emd_D);
    [MSE, PSNR, ER, BPP] = caulate(row, col, inter_I, recover_I, num_emd);

    disp(['image:' I_name])
    disp(['Ƕ�������: ' num2str(num_emd) ' bits'] )
    disp(['Ƕ����: ' num2str(BPP) ' bpp'])
    disp(['������: ' num2str(ER*100) '%'])

end