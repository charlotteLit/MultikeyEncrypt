I = imread('./测试图像/Lena.tiff');
origin_I = double(I);
[row,col] = size(origin_I);

double_I = [origin_I, origin_I];
% figure(1);
% imshow(double_I,[]);
% title('拼接图像');

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
d = 0;
encrypt_I = double_I;

% 第一行第一列灰集
% 加密灰集，选定密钥1
key = squeeze(rand_seq(1,:,:));
encrypt_I(1:bs,:) = bitxor(encrypt_I(1:bs,:), key(1:bs,:));
encrypt_I(bs+1:510,1:bs) = bitxor(encrypt_I(bs+1:510,1:bs),key(bs+1:510,1:bs));
encrypt_I(bs+1:510,511:515) = bitxor(encrypt_I(bs+1:510,511:515),key(bs+1:510,511:515));
encrypt_I(bs+1:510,end-2:end) = bitxor(encrypt_I(bs+1:510,end-2:end),key(bs+1:510,end-2:end));
encrypt_I(510:512,:) = bitxor(encrypt_I(510:512,:),key(510:512,:));


figure(1);
imshow(encrypt_I,[]);
title('加密图像1');

emd_data = zeros();
emd = 0;

% 加密白集
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
        % 对第一幅图像的相应块进行加密
        encrypt_I(x_range, y_range) = bitxor(encrypt_I(x_range, y_range), k1);
        % 对第二幅图像的相应块进行加密
        encrypt_I(x_range, y1_range) = bitxor(encrypt_I(x_range, y1_range), k2);
    end
end

figure(2);
imshow(encrypt_I,[]);
title('加密图像2');

% 数据提取 & 复原图像


