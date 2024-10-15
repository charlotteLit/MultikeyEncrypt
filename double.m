I = imread('./����ͼ��/Lena.tiff');
origin_I = double(I);
[row,col] = size(origin_I);

double_I = [origin_I, origin_I];
% figure(1);
% imshow(double_I,[]);
% title('ƴ��ͼ��');

% �ֿ��С
bs = 3;

% �ֿ���
m = floor(row/bs);
n = floor(col/bs);

% Ƕ����������
num_D = 1000000;
rng(0);
D = round(rand(1,num_D)*1);

% ����2^b��������У���СΪ����ͼƬˮƽ����
% ÿ��Ƕ��5 bit����
b = 5;
rand_num = 2^b;
rng(b);
rand_seq = floor(rand(rand_num, row, col*2) * 255);

% ͼ����� & ����Ƕ��
% ������Ϣ���
d = 0;
encrypt_I = double_I;

% ��һ�е�һ�лҼ�
% ���ܻҼ���ѡ����Կ1
key = squeeze(rand_seq(1,:,:));
encrypt_I(1:bs,:) = bitxor(encrypt_I(1:bs,:), key(1:bs,:));
encrypt_I(bs+1:510,1:bs) = bitxor(encrypt_I(bs+1:510,1:bs),key(bs+1:510,1:bs));
encrypt_I(bs+1:510,511:515) = bitxor(encrypt_I(bs+1:510,511:515),key(bs+1:510,511:515));
encrypt_I(bs+1:510,end-2:end) = bitxor(encrypt_I(bs+1:510,end-2:end),key(bs+1:510,end-2:end));
encrypt_I(510:512,:) = bitxor(encrypt_I(510:512,:),key(510:512,:));


figure(1);
imshow(encrypt_I,[]);
title('����ͼ��1');

emd_data = zeros();
emd = 0;

% ���ܰ׼�
for i = 2:m
    begin_x = (i-1) * bs;
    x_range = begin_x+1:begin_x+bs;
    
    for j = 2:n
        begin_y = (j-1) * bs;
        begin_y1 = begin_y + col;
        y_range = begin_y+1:begin_y+bs;
        y1_range = begin_y1+1:begin_y1+bs;

        ed = bin2dec([0, 0, 0, D(d+1:d+5)]) + 1;
        d = d + 5;
        emd_data(emd+1:emd+5) = D(d+1:d+5);
        emd = emd + 5;

        key = rand_seq(ed, :, :);
        key = reshape(key, row, col*2);
        k1 = key(x_range, y_range);
        k2 = key(x_range, y1_range);
        % �Ե�һ��ͼ�����Ӧ����м���
        encrypt_I(x_range, y_range) = bitxor(encrypt_I(x_range, y_range), k1);
        % �Եڶ���ͼ�����Ӧ����м���
        encrypt_I(x_range, y1_range) = bitxor(encrypt_I(x_range, y1_range), k2);
    end
end

figure(2);
imshow(encrypt_I,[]);
title('����ͼ��2');

% ������ȡ & ��ԭͼ��


