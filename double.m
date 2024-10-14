I = imread('./����ͼ��/Lena.tiff');
origin_I = double(I);
[row,col] = size(origin_I);

double_I = [origin_I, origin_I];
% figure(1);
% imshow(double_I,[]);title('ƴ��ͼ��');

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
d = 1;
encrypt_I = double_I;
% ��һ�е�һ�лҼ�
% ���ܻҼ���ѡ����Կ1
key = rand_seq(1,:,:);
encrypt_I(1:bs,:) = bitxor(encrypt_I(1:bs,:), key(1:bs,:));
encrypt_I(bs+1:end,1:bs) = bitxor(encrypt_I(bs+1:end,1:bs),key(bs+1:end,1:bs));

% ���ܰ׼�
for i=2:m
    for j=2:n
        begin_x = (i-1)*bs+1;
        begin_y = (j-1)*bs+1;
        begin_y1 = begin_y + col;
        ed = bin2dec([0,0,0,D(d:d+5)]);
        d = d+5;
        key = rand_seq(ed,:,:);
        % ���ܵ�һ��ͼ���Ӧ��
        encrypt_I(begin_x:begin_x+bs,begin_y:begin_y+bs) = bitxor(encrypt_I(begin_x:begin_x+bs,begin_y:begin_y+bs), ...
        key(begin_x:begin_x+bs,begin_y:begin_y+bs));
        % ���ܵڶ���ͼ���Ӧ��
        encrypt_I(begin_x:begin_x+bs,begin_y1:begin_y1+bs) = bitxor(encrypt_I(begin_x:begin_x+bs,begin_y1:begin_y1+bs), ...
        key(begin_x:begin_x+bs,begin_y1:begin_y1+bs));
    end
end

% ������ȡ & ��ԭͼ��


