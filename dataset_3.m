clc
clear

% �ֿ��С
bs = 3;

% Ƕ����������
num_D = 1000000;
rng(0);
D = round(rand(1,num_D)*1);

% ÿ��Ƕ��5 bit����
b = 5;

%% ͼ�����ݼ���Ϣ(BOSSbase_1.01),��ʽ:PGM,����:10000��
I_file_path = 'D:\ImageDatabase\BOSSbase_1.01\'; %����ͼ�����ݼ��ļ���·��
I_path_list = dir(strcat(I_file_path,'*.pgm'));  %��ȡ���ļ���������pgm��ʽ��ͼ��
img_num = length(I_path_list); %��ȡͼ��������
%% ��¼ÿ��ͼ��ĵ������Ϣ
bpp_BOSSbase = zeros(1,img_num); %��¼ÿ��ͼ���Ƕ����
er_BOSSbase = zeros(1,img_num);%��¼ÿ��ͼ���������ظ���

for i=1:img_num
    %-------------------------------��ȡͼ��-------------------------------%
    I_name = I_path_list(i).name; %ͼ����
    I = imread(strcat(I_file_path,I_name));%��ȡͼ��
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
        disp('��ȡ������Ƕ��������ȫ��ͬ��')
        disp(['Ƕ�������: ' num2str(num_emD) ' bits'] )
        disp(['Ƕ����: ' num2str(bpp) ' bpp'])
    else
        disp('Warning��������ȡ����')
    end
    fprintf(['�� ',num2str(i),' ��ͼ��-------- OK','\n\n']);
end
