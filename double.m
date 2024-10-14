I = imread('./测试图像/Lena.tiff');
origin_I = double(I);
[row,col] = size(origin_I);

double_I = [origin_I, origin_I];
% figure(1);
% imshow(double_I,[]);title('拼接图像');

% 分块大小
bs = 3;

% 分块数
m = floor(row/bs);
n = floor(col/bs);

% 嵌入数据生成
num_D = 1000000;
rng(0);
D = round(rand(1,num_D)*1);

% 生成2^b个随机序列，大小为两幅图片水平排列
% 每块嵌入5 bit数据
b = 5;
rand_num = 2^b;
rng(b);
rand_seq = floor(rand(rand_num, row, col*2) * 255);

% 图像加密 & 数据嵌入
% 隐藏信息标记
d = 1;
encrypt_I = double_I;
% 第一行第一列灰集
% 加密灰集，选定密钥1
key = rand_seq(1,:,:);
encrypt_I(1:bs,:) = bitxor(encrypt_I(1:bs,:), key(1:bs,:));
encrypt_I(bs+1:end,1:bs) = bitxor(encrypt_I(bs+1:end,1:bs),key(bs+1:end,1:bs));

% 加密白集
for i=2:m
    for j=2:n
        begin_x = (i-1)*bs+1;
        begin_y = (j-1)*bs+1;
        begin_y1 = begin_y + col;
        ed = bin2dec([0,0,0,D(d:d+5)]);
        d = d+5;
        key = rand_seq(ed,:,:);
        % 加密第一幅图像对应块
        encrypt_I(begin_x:begin_x+bs,begin_y:begin_y+bs) = bitxor(encrypt_I(begin_x:begin_x+bs,begin_y:begin_y+bs), ...
        key(begin_x:begin_x+bs,begin_y:begin_y+bs));
        % 加密第二幅图像对应块
        encrypt_I(begin_x:begin_x+bs,begin_y1:begin_y1+bs) = bitxor(encrypt_I(begin_x:begin_x+bs,begin_y1:begin_y1+bs), ...
        key(begin_x:begin_x+bs,begin_y1:begin_y1+bs));
    end
end

% 数据提取 & 复原图像


